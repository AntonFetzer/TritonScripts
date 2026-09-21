"""Validate and inspect the HUS production water response matrix and exploratory fits."""
from pathlib import Path
import csv, json, re, sys, hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from Dependencies.PoolDoseModules import poolDoseModules
from Read.ReadDose import readDoseModules
from Dependencies.ScanMetrics import cellBounds, fallingCrossing
from Read.ReadMCC import readScan
from Plotting.RadExHUSScans import DEPTHS_MM

BASE=Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS")
RUN=BASE/"PDD-Spectrum-Production"

def _csv(path, rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

def main():
    out=RUN/"Analysis"; out.mkdir(exist_ok=True)
    # Historical raw outputs remain unchanged. Hashes prevent applying this
    # correction to replacement CURRENT-normalized results or mixed datasets.
    normalization=json.loads((RUN/"normalization.json").read_text())
    factor=float(normalization["raw_to_current_factor"])
    if normalization["output_convention"] != "CURRENT" or factor not in (1.0,4.0):
        raise ValueError("Unsupported normalization provenance")
    raw_files={str(p.relative_to(RUN)):p for p in (RUN/"Res").rglob("*.csv")}
    if set(raw_files) != set(normalization["raw_sha256"]):
        raise ValueError("Raw dataset differs from normalization provenance")
    for name,p in raw_files.items():
        if hashlib.sha256(p.read_bytes()).hexdigest() != normalization["raw_sha256"][name]:
            raise ValueError(f"Raw normalization hash mismatch: {p}")
    pilot_path=out/"pilot_reference.csv"
    if hashlib.sha256(pilot_path.read_bytes()).hexdigest() != normalization["pilot_reference_sha256"]:
        raise ValueError("Pilot normalization hash mismatch")
    pilot_factor=float(normalization["pilot_to_current_factor"])
    manifest=list(csv.DictReader((RUN/"manifest.csv").open()))
    energies=sorted({float(r["energy_MeV"]) for r in manifest})
    K=[]; S=[]; NZ=[]; B=[]; rows=[]; seeds=set(); expected=set()
    for e in energies:
        files=[]
        for r in [r for r in manifest if float(r["energy_MeV"])==e]:
            prefix=RUN/r["output_prefix"]
            matches=list(prefix.parent.glob(prefix.name+"_*.csv"))
            if len(matches)!=1: raise ValueError(f"Missing/duplicate replica {r}")
            p=matches[0]; expected.add(p); files.append(p)
            pair=tuple(p.stem.split("_")[-2:])
            if pair in seeds: raise ValueError(f"Repeated seed pair {pair}")
            seeds.add(pair)
            d=readDoseModules(p,"doseLayer-")
            if set(d)!={f"doseLayer-{i}" for i in range(1,35)}: raise ValueError("Modules")
            if any(v["entries"]!=int(r["histories"]) or v["unit"]!="MeV/g" for v in d.values()):
                raise ValueError(f"Entries/units: {p}")
        if len(files)!=5: raise ValueError("Replica count")
        pooled=poolDoseModules(files)
        vals=[pooled[f"doseLayer-{i}"] for i in range(1,35)]
        for value in vals:
            value["dose"] *= factor
            value["error"] *= factor
        K.append([v["dose"] for v in vals]); S.append([v["error"] for v in vals])
        NZ.append([v["non-zeros"] for v in vals]); B.append([v["birge_ratio"] for v in vals])
        for z,v in zip(DEPTHS_MM,vals):
            rows.append(dict(energy_MeV=e,depth_mm=z,**v))
    if expected!=set((RUN/"Res").rglob("*.csv")): raise ValueError("Unexpected raw results")
    K=np.array(K).T; S=np.array(S).T; NZ=np.array(NZ).T; B=np.array(B).T
    energies=np.array(energies)
    rel=np.divide(100*S,K,out=np.full_like(K,np.nan),where=K>0)
    pdd=100*K/K.max(axis=0)
    _csv(out/"pooled_kernels.csv",rows)
    np.savez(out/"response_matrix.npz",energy_MeV=energies,depth_mm=DEPTHS_MM,
             dose=K,error=S,entries_per_energy=np.full(len(energies),5000000),
             nonzero=NZ,birge=B,unit="MeV/g at incident CURRENT fluence 1 cm^-2",
             raw_to_current_factor=factor)
    for name,a in (("response_matrix.csv",K),("response_matrix_error.csv",S)):
        _csv(out/name,[dict(depth_mm=z,**{f"E{e:.2f}_MeV":a[i,j] for j,e in enumerate(energies)}) for i,z in enumerate(DEPTHS_MM)])
    metrics=[]
    for j,e in enumerate(energies):
        active=pdd[:,j]>=1
        metrics.append(dict(energy_MeV=e,dmax_mm=DEPTHS_MM[np.argmax(K[:,j])],
            R50_mm=fallingCrossing(DEPTHS_MM,pdd[:,j],50),
            surface_percent=pdd[0,j],max_relative_error_percent=float(np.nanmax(rel[:,j])),
            max_error_above_1percent_peak=float(np.nanmax(rel[active,j])),
            error_60mm_percent=rel[-1,j],max_birge=B[:,j].max()))
    _csv(out/"kernel_metrics.csv",metrics)
    # Retain only the derived reference table after removal of the pilot runs.
    reference=list(csv.DictReader((out/"pilot_reference.csv").open()))
    pilot=[]
    for r in reference:
        e=float(r["energy_MeV"]); z=float(r["depth_mm"])
        j=np.where(energies==e)[0][0]; i=np.where(DEPTHS_MM==z)[0][0]
        sigma=np.hypot(S[i,j],pilot_factor*float(r["error_MeV_per_g_per_cm2"]))
        pilot.append(dict(energy_MeV=e,depth_mm=z,
            z_score=(K[i,j]-pilot_factor*float(r["dose_MeV_per_g_per_cm2"]))/sigma))
    _csv(out/"pilot_comparison.csv",pilot)
    scan=readScan(next(BASE.glob("*eHDTSE_PDD.mcc")),"PDD")
    depths,measurement=scan["position"],scan["value"]
    if not np.array_equal(depths,DEPTHS_MM): raise ValueError("Depth mismatch")
    measured=measurement/measurement.max()
    # Equal absolute PDD weights: no measurement error model is supplied.
    # Absolute kernels share ONE global scale; no per-energy peak normalization.
    A=K/K.max()
    smooth=np.diff(np.eye(len(energies)),n=2,axis=0)
    fits=[]; solutions=[]
    for limit in (30.,40.):
        mask=depths<=limit
        for lam in (0.,1e-4,1e-3,1e-2,1e-1,1.):
            x,_=nnls(np.vstack((A[mask],np.sqrt(lam)*smooth)),
                     np.r_[measured[mask],np.zeros(len(smooth))],maxiter=10000)
            pred=A@x; weights=x/x.sum()
            fits.append(dict(depth_limit_mm=limit,regularization=lam,
                rmse_fit_pp=100*np.sqrt(np.mean((pred[mask]-measured[mask])**2)),
                max_abs_residual_fit_pp=100*np.max(np.abs(pred[mask]-measured[mask])),
                mean_energy_MeV=float(energies@weights),
                peak_energy_MeV=float(energies[np.argmax(weights)]),
                fraction_below_2MeV=float(weights[energies<2].sum()),
                fraction_above_7MeV=float(weights[energies>7].sum()),
                predicted_surface_percent=100*pred[0],predicted_60mm_percent=100*pred[-1]))
            solutions.append((limit,lam,weights,pred,x))
    _csv(out/"exploratory_fit_metrics.csv",fits)
    _csv(out/"exploratory_spectra.csv",[dict(depth_limit_mm=l,regularization=a,energy_MeV=e,fluence_fraction=w[j])
        for l,a,w,p,x in solutions for j,e in enumerate(energies)])
    _csv(out/"exploratory_predictions.csv",[dict(depth_limit_mm=l,regularization=a,depth_mm=z,
        predicted_percent=100*p[i],measured_percent=100*measured[i],residual_pp=100*(p[i]-measured[i]))
        for l,a,w,p,x in solutions for i,z in enumerate(depths)])
    flags=[r for r in rows if r["inconsistent"]]
    summary=dict(normalization=normalization["output_convention"],
        raw_to_current_factor=factor, pilot_to_current_factor=pilot_factor, files=len(expected),energies=len(energies),histories=160000000,
        unique_seed_pairs=len(seeds),max_birge=float(B.max()),birge_flags=flags,
        max_pilot_abs_z=max(abs(r["z_score"]) for r in pilot),
        max_relative_error_percent=float(np.nanmax(rel)),
        max_error_above_1percent_peak=max(r["max_error_above_1percent_peak"] for r in metrics),
        measured_surface_percent=100*measured[0],measured_60mm_percent=100*measured[-1])
    (out/"summary.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    print(json.dumps(fits,indent=2))
    # One multipage report; relative-error axes and color scale are linear.
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out/"PDD-production-analysis.pdf") as pdf:
        fig,axs=plt.subplots(2,1,figsize=(8,9))
        for e in (0.25,1,2,3,4,5,5.5,6,6.5,7,8):
            j=np.where(energies==e)[0][0]
            axs[0].plot(depths,pdd[:,j],label=f"{e:g} MeV",lw=1)
        axs[0].plot(depths,100*measured,"ko",ms=3,label="HUS measured")
        axs[0].set(xlabel="Depth [mm]",ylabel="Dose / kernel maximum [%]",title="Production water kernels")
        axs[0].legend(ncol=3,fontsize=8);axs[0].grid(alpha=.3)
        im=axs[1].imshow(rel.T,aspect="auto",origin="lower",extent=(-.75,62.5,.125,8.125),vmin=0)
        axs[1].set(xlabel="Depth [mm] (irregular samples shown in index order)",ylabel="Energy [MeV]",title="Relative statistical error [%], linear scale")
        # pcolormesh respects the irregular measured-depth cells.
        axs[1].clear()
        im=axs[1].pcolormesh(cellBounds(depths),np.arange(.125,8.126,.25),rel.T,vmin=0,shading="flat")
        axs[1].set(xlabel="Depth [mm]",ylabel="Energy [MeV]",title="Relative statistical error [%], linear scale")
        fig.colorbar(im,ax=axs[1],label="Relative error [%]")
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        fig,axs=plt.subplots(2,1,figsize=(8,8))
        for e in (0.25,1,2,4,5,5.5,6,8):
            j=np.where(energies==e)[0][0]
            axs[0].plot(depths,rel[:,j],"o-",ms=2,label=f"{e:g} MeV")
            axs[1].plot(depths,np.where(pdd[:,j]>=1,rel[:,j],np.nan),"o-",ms=2,label=f"{e:g} MeV")
        for ax in axs:
            ax.set(xlabel="Depth [mm]",ylabel="Relative statistical error [%]",ylim=(0,None))
            ax.axhline(1,color="k",ls="--",lw=1);ax.grid(alpha=.3);ax.legend(ncol=4,fontsize=8)
        axs[0].set_title("All depths: high errors occur in weak low-energy tails")
        axs[1].set_title("Only cells with dose >= 1% of that kernel's maximum")
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        fig,axs=plt.subplots(3,1,figsize=(8,10))
        axs[0].plot(depths,100*measured,"ko",ms=3,label="HUS measured")
        for limit,lam,w,p,x in solutions:
            if lam in (0.,.01,.1) and limit==40.:
                label=f"lambda={lam:g}, fit <=40 mm"
                axs[0].plot(depths,100*p,label=label)
                axs[1].plot(depths,100*(p-measured),label=label)
                axs[2].plot(energies,100*w,"o-",ms=3,label=label)
        axs[0].set(ylabel="PDD [%]",title="Exploratory electron-only fits; not a commissioned spectrum")
        axs[1].set(xlabel="Depth [mm]",ylabel="Fit - measurement [pp]")
        axs[2].set(xlabel="Electron energy [MeV]",ylabel="Fluence fraction per 0.25 MeV [%]")
        for ax in axs: ax.grid(alpha=.3);ax.legend(fontsize=8)
        fig.text(.5,.01,"Equal absolute PDD weights; measurement uncertainties unavailable.\nSpectral variation across regularization is sensitivity, not a confidence interval.",ha="center",fontsize=8)
        fig.tight_layout(rect=(0,.045,1,1));pdf.savefig(fig);plt.close(fig)
        fig,axs=plt.subplots(2,1,figsize=(8,8))
        for e in (5.,5.5,6.):
            pp=[r for r in pilot if r["energy_MeV"]==e]
            axs[0].plot(depths,[r["z_score"] for r in pp],"o-",ms=3,label=f"{e:g} MeV")
        axs[0].axhline(0,color="k",lw=.7);axs[0].set(xlabel="Depth [mm]",ylabel="Production - pilot [combined sigma]",title="Independent pilot comparison")
        for limit,lam,w,p,x in solutions:
            if lam==.01:
                axs[1].plot(energies,100*w,"o-",label=f"Fit through {limit:g} mm")
        axs[1].set(xlabel="Electron energy [MeV]",ylabel="Fluence fraction [%]",title="Fit-window sensitivity at lambda=0.01")
        for ax in axs:ax.grid(alpha=.3);ax.legend()
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
if __name__=="__main__":
    main()
