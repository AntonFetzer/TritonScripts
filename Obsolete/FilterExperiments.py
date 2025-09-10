import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.odr as odr
import pandas as pd
from GRAS.Dependencies.TotalDose import totalkRadGras


def simpleFilter(data, window_size=5):
    dose = pd.Series(data[0])
    error = pd.Series(data[1])
    dose_windows = dose.rolling(window_size, center=True)
    error_windows = error.rolling(window_size, center=True)

    errors = []
    for window in error_windows:
        squares = window ** 2
        new_error = np.sqrt(squares.sum()) / window_size
        errors.append(new_error)

    mov_avg = dose_windows.mean()
    final_avg = data[0][:int(window_size / 2)].tolist() + mov_avg.tolist()[int(window_size / 2):-int(window_size / 2)] + \
                data[0][-int(window_size / 2):].tolist()
    final_error = data[1][:int(window_size / 2)].tolist() + errors[int(window_size / 2):-int(window_size / 2)] + data[
                                                                                                                     1][
                                                                                                                 -int(
                                                                                                                     window_size / 2):].tolist()

    # for i in range(len(data) - 4):
    #     Res[0][i + 2] = (data[0][i] + data[0][i + 1] + data[0][i + 2] + data[0][i + 3] + data[0][i + 4]) / 5
    #     # print(Res[i], data[i])
    #     Res[1][i + 2] = np.sqrt(data[1][i]**2 + data[1][i + 1]**2 + data[1][i + 2]**2 + data[1][i + 3]**2 + data[1][i + 4]**2) / 5

    return np.array([final_avg, final_error])


if __name__ == "__main__":
    path = "/l/triton_work/2LayerStackedCurves/PE-Pb-32/Res/"

    # Data = readGrasCsv(path)
    Data = totalkRadGras(path, "Elec")

    x = np.linspace(0, 100, num=101, dtype=int)
    #plt.plot(x, Data[0], 'C0.')  # , label="Raw")
    plt.errorbar(x, Data[0], Data[1], fmt='C0.', capsize=2, label="Raw")
    plt.plot(np.argmin(Data[0]), min(Data[0]), "C0X",
             label="Raw Min=%5.2f at %d %%" % (min(Data[0]), np.argmin(Data[0])))

    WindowSize = 10
    SimpleFiltered = simpleFilter(Data, WindowSize)

    # for i in range(5):
    #    SimpleFiltered = simpleFilter(SimpleFiltered)
    plt.errorbar(x, SimpleFiltered[0], SimpleFiltered[1], fmt="C1", capsize=2)  # , label="SimpleFiltered")
    MinX = np.nanargmin(SimpleFiltered[0])
    MinY = np.nanmin(SimpleFiltered[0])
    plt.plot(MinX, MinY, "C1X", label="SimpleFiltered Min=%5.2f at %d %%" % (MinY, MinX))

    '''
    savgol = Data[0]
    for i in range(1):
        savgol = signal.savgol_filter(savgol, 10, 4)
    plt.plot(x, savgol, "C2", label="savgol")
    plt.plot(np.argmin(savgol), min(savgol), "C2X")


    WindowSize = 10
    Fitsubset = Data[0][MinX-WindowSize:MinX+WindowSize]
    xSubset = x[MinX-WindowSize:MinX+WindowSize]
    Fit = np.poly1d(np.polyfit(xSubset, Fitsubset, 2))
    Fit = Fit(x)
    plt.plot(xSubset, Fitsubset, "C3o", label="Fitsubset")
    plt.plot(x, Fit, "C3", label="Fit")
    plt.plot(np.argmin(Fit), min(Fit), "C3X")
    '''

    WindowStart = MinX - WindowSize
    WindowEnd = MinX + WindowSize

    if WindowSize > MinX:
        WindowStart = 0
        WindowEnd = 2 * WindowSize
    elif WindowSize > len(Data[0]) - MinX:
        WindowEnd = len(Data[0])
        WindowStart = len(Data[0]) - 2 * WindowSize

    #FitData = [SimpleFiltered[0][WindowStart:WindowEnd], SimpleFiltered[1][WindowStart:WindowEnd]]
    FitData = [Data[0][WindowStart:WindowEnd], Data[1][WindowStart:WindowEnd]]
    xSet = x[WindowStart:WindowEnd]
    plt.errorbar(xSet, FitData[0], FitData[1], fmt='C3.', capsize=2)  # , label="FitData")


    def func(x, a, b, c):
        return a * x * x + b * x + c


    # fitting with curve with only mean square optimsation
    parametererLimits = ([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf])
    # popt, pcov = curve_fit(func, xSet, FitData[0])
    popt, pcov = curve_fit(func, xSet, FitData[0], sigma=FitData[1], bounds=parametererLimits, absolute_sigma=True)

    FitMeanSquare = func(x, *popt)
    plt.plot(x, FitMeanSquare, 'C3')  # , label="Fit")
    SimpleMinX = np.argmin(FitMeanSquare)
    SimpleMinY = min(FitMeanSquare)

    perr = np.sqrt(np.diag(pcov))
    print(perr)

    Cfit = popt[2]
    CErrfit = perr[2]
    Bfit = popt[1]
    BErrfit = perr[1]
    Afit = popt[0]
    AErrfit = perr[0]
    FitMinX = -Bfit / (2 * Afit)
    FitMinY = Cfit - (Bfit ** 2) / (4 * Afit)

    SimpleYError = np.sqrt( (SimpleMinX**2 * AErrfit)**2 + (SimpleMinX * BErrfit)**2 + CErrfit**2)
    plt.errorbar(SimpleMinX, SimpleMinY, color="C5", yerr=SimpleYError, fmt="X",
             label="Fit Min=%5.2f (%.2f) at %d (%d) %%" % (SimpleMinY, SimpleMinY, SimpleMinY, 0))

    FitMinXErr = np.sqrt((AErrfit * Bfit / (2 * Afit ** 2)) ** 2 + (BErrfit / (2 * Afit)) ** 2)
    FitMinYErr = np.sqrt(CErrfit ** 2 + (Bfit * BErrfit / (2 * Afit)) ** 2 + (AErrfit * Bfit ** 2 / (4 * Afit ** 2)) ** 2)
    plt.errorbar(FitMinX, FitMinY, color="C3", xerr=FitMinXErr, yerr=FitMinYErr, fmt="X",
             label="Fit Min=%5.2f (%.2f) at %d (%d) %%" % (FitMinY, FitMinYErr, FitMinX, FitMinXErr))

    # fitting with orthogonal distance regression (also error in x)
    parabola = odr.Model(func)
    mydata = odr.RealData(xSet, FitData[0], sy=FitData[1])
    myodr = odr.ODR(mydata, odr.quadratic)  # , beta0=[1, 1, MinY])
    myoutput = myodr.run()
    #myoutput.pprint()

    print("myoutput.sum_square: ", myoutput.sum_square)
    print("myoutput.sum_square_delta: ", myoutput.sum_square_delta)

    FitOdr = func(x, *myoutput.beta)
    plt.plot(x, FitOdr, 'C2')  # , label="Fit ODR")
    COdr = myoutput.beta[2]
    CErrOdr = myoutput.sd_beta[2]
    BOdr = myoutput.beta[1]
    BErrOdr = myoutput.sd_beta[1]
    AOdr = myoutput.beta[0]
    AErrOdr = myoutput.sd_beta[0]
    OdrMinX = -BOdr / (2 * AOdr)
    OdrMinY = COdr - (BOdr ** 2) / (4 * AOdr)

    FitMinXErrOdr = np.sqrt((AErrOdr * BOdr / (2 * AOdr ** 2)) ** 2 + (BErrOdr / (2 * AOdr)) ** 2)
    FitMinYErrOdr = np.sqrt(CErrOdr ** 2 + (BOdr * BErrOdr / (2 * AOdr)) ** 2 + (AErrOdr * BOdr ** 2 / (4 * AOdr ** 2)) ** 2)

    plt.errorbar(OdrMinX, OdrMinY, color="C2", xerr=FitMinXErrOdr, yerr=FitMinYErrOdr, fmt="X",
                 label="FitOdr Min=%5.2f (%.2f) at %d (%d) %%" % (OdrMinY, FitMinYErrOdr, OdrMinX, FitMinXErrOdr))

    plt.grid(which='both')
    plt.legend()
    # plt.yscale("log")
    #plt.ylim(0, 0.02)
    plt.show()
