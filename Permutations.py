Materials = ["Polyethylene", "Aluminium", "Tungsten", "Lead"]

for i1 in range(4):
    for i2 in range(4):
        for i3 in range(4):
            for i4 in range(4):
                for i5 in range(4):
                    print(Materials[i5], Materials[i4], Materials[i3], Materials[i2], Materials[i1])

# print(Shields[999])
'''
Prot = np.loadtxt("/home/anton/Desktop/triton_work/Permutations/5Layer1-25gcm2/csv/5layer2e7proton.txt")
Elec = np.loadtxt("/home/anton/Desktop/triton_work/Permutations/5Layer1-25gcm2/csv/5layer2e9electron.txt")


NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 2e9 / 1024
Npart_Prot = 2e7 / 1024

for i in range(1297):
    print(MeVtokRad_2D(Prot[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot), MeVtokRad_2D(Elec[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec))
'''