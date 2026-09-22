#!/usr/bin/env python3
import csv, math, os, re, statistics
from pathlib import Path

ROOT=Path("pressure_history_extract")
OUT=Path("analysis")
OUT.mkdir(exist_ok=True)

SOLVER_ALIASES = {
    "RMT":"RMT",
    "RMT_FINALVOL2":"RMT",
    "MG3V":"MG3V",
    "MG3V_3LEVEL":"MG3V",
    "MG2V":"MG2V",
    "MG2V_3LEVEL":"MG2V",
    "MG2W":"MG2W",
    "MG2W_3LEVEL":"MG2W",
    "SG":"SG",
    "SG_RBGS":"SG",
}
WINDOWS=[0.001,0.005,0.01,0.05,0.1,0.25,0.5,1.0,2.0]

def fnum(x):
    try:
        return float(x)
    except Exception:
        return math.nan

def infer_solver(path, header):
    text=str(path).upper()
    for key in ["RMT_FINALVOL2","RMT","MG3V","MG2W","MG2V","SG_RBGS","SG"]:
        if key in text:
            return SOLVER_ALIASES[key]
    # Header convention itself cannot distinguish MG families, so use file/path first.
    return "UNKNOWN"

# Campaign rows at the requested physical condition: intermediate stiffness, max mesh, max core.
cost={}
campaign_specs=[
    ("RMT","RMT_FINALVOL2_campaign_summary.csv","total_rmt_cycles"),
    ("MG3V","MG3V_campaign_summary.csv","total_vcycles"),
    ("MG2V","MG2V_campaign_summary.csv","total_vcycles"),
    ("MG2W","MG2W_campaign_summary.csv","total_wcycles"),
    ("SG","SG_RBGS_campaign_summary.csv","total_rbgs_sweeps"),
]
condition_rows=[]
for solver, fn, totalcol in campaign_specs:
    p=Path(fn)
    if not p.exists():
        continue
    with p.open(newline="") as fh:
        rows=list(csv.DictReader(fh))
    cand=[]
    for r in rows:
        mesh=r.get("mesh_name","")
        th=r.get("threads","")
        pc=r.get("physical_case","")
        stiff=r.get("stiffness_class","")
        # Exact requested condition if present.
        if mesh=="M4_864x288" and th=="16" and (
            "S2_INTERMEDIATE" in pc or "INTERMEDIATE" in pc or stiff=="intermediate"
        ):
            cand.append(r)
    if not cand:
        # Fallback: max mesh/core and intermediate-like reaction case.
        for r in rows:
            if r.get("mesh_name","")=="M4_864x288" and r.get("threads","")=="16" and "INTERMEDIATE" in (r.get("case_id","")+r.get("physical_case","")).upper():
                cand.append(r)
    if cand:
        r=cand[0]
        pt=fnum(r.get("pressure_time_s","nan"))
        nt=fnum(r.get(totalcol,"nan"))
        if math.isfinite(pt) and math.isfinite(nt) and nt>0:
            cost[solver]=pt/nt
        condition_rows.append({
            "solver":solver,
            "case_id":r.get("case_id",""),
            "physical_case":r.get("physical_case",""),
            "mesh_name":r.get("mesh_name",""),
            "threads":r.get("threads",""),
            "pressure_time_s":pt,
            "total_work_units":nt,
            "seconds_per_work_unit":cost.get(solver,math.nan),
            "total_wall_s":fnum(r.get("total_wall_s","nan")),
        })

# Inventory and per-history analysis.
files=sorted([p for p in ROOT.rglob("*.csv") if "pressure_step_history" in p.name.lower() or "pressure" in p.name.lower()])
inventory=[]
all_rows=[]
selected_files=[]

for p in files:
    try:
        with p.open(newline="", errors="replace") as fh:
            rd=csv.DictReader(fh)
            header=rd.fieldnames or []
            rows=list(rd)
    except Exception as e:
        inventory.append((str(p),0,"READ_ERROR",str(e)))
        continue
    if not rows:
        inventory.append((str(p),0,"EMPTY",""))
        continue
    solver=infer_solver(p,header)
    inventory.append((str(p),len(rows),solver,",".join(header)))
    # Require timestep history fields.
    tcol=next((x for x in ["time_s","time"] if x in header),None)
    stepcol=next((x for x in ["step"] if x in header),None)
    cyc_col=next((x for x in ["cycles","rmt_cycles","vcycles","wcycles","rbgs_sweeps","sweeps"] if x in header),None)
    if cyc_col is None:
        # Production files preserve a legacy solver-specific fourth header.
        if len(header)>=4:
            cyc_col=header[3]
    init_col=next((x for x in ["initial_residual","initialResidual"] if x in header),None)
    final_col=next((x for x in ["final_residual","finalResidual"] if x in header),None)
    rel_col=next((x for x in ["final_relative_residual","relative_residual","final_rel_residual"] if x in header),None)
    rho_col=next((x for x in ["convergence_factor","rho","residual_ratio"] if x in header),None)
    if tcol is None or cyc_col is None:
        continue
    parsed=[]
    for r in rows:
        t=fnum(r.get(tcol,"nan")); c=fnum(r.get(cyc_col,"nan"))
        if not (math.isfinite(t) and math.isfinite(c)):
            continue
        parsed.append({
            "t":t,
            "step":fnum(r.get(stepcol,"nan")) if stepcol else math.nan,
            "c":c,
            "init":fnum(r.get(init_col,"nan")) if init_col else math.nan,
            "final":fnum(r.get(final_col,"nan")) if final_col else math.nan,
            "rel":fnum(r.get(rel_col,"nan")) if rel_col else math.nan,
            "rho":fnum(r.get(rho_col,"nan")) if rho_col else math.nan,
        })
    if not parsed:
        continue
    selected_files.append((solver,p,len(parsed),cyc_col))
    sec=cost.get(solver,math.nan)
    for w in WINDOWS:
        rr=[x for x in parsed if x["t"]<=w+1e-15]
        if not rr:
            continue
        cs=[x["c"] for x in rr]
        inits=[x["init"] for x in rr if math.isfinite(x["init"])]
        rels=[x["rel"] for x in rr if math.isfinite(x["rel"])]
        all_rows.append({
            "solver":solver,
            "source_file":str(p),
            "window_end_s":w,
            "solves":len(rr),
            "cumulative_work_units":sum(cs),
            "mean_work_units_per_solve":sum(cs)/len(cs),
            "zero_work_fraction":sum(1 for x in cs if x==0)/len(cs),
            "one_or_less_fraction":sum(1 for x in cs if x<=1)/len(cs),
            "max_work_units":max(cs),
            "initial_residual_median":statistics.median(inits) if inits else math.nan,
            "initial_residual_max":max(inits) if inits else math.nan,
            "final_rel_residual_max":max(rels) if rels else math.nan,
            "sec_per_work_unit_from_campaign":sec,
            "estimated_cumulative_pressure_time_s":sum(cs)*sec if math.isfinite(sec) else math.nan,
        })

# Write inventory.
with (OUT/"pressure_history_inventory.txt").open("w") as fh:
    fh.write("Extracted pressure-history inventory\n")
    fh.write("===================================\n")
    for row in inventory:
        fh.write(" | ".join(map(str,row))+"\n")
    fh.write("\nSelected histories\n")
    for solver,p,n,cc in selected_files:
        fh.write(f"{solver} | {p} | rows={n} | work_col={cc}\n")
    fh.write("\nCampaign cost anchors\n")
    for r in condition_rows:
        fh.write(str(r)+"\n")

# Detailed window summary.
fields=[
    "solver","source_file","window_end_s","solves","cumulative_work_units",
    "mean_work_units_per_solve","zero_work_fraction","one_or_less_fraction",
    "max_work_units","initial_residual_median","initial_residual_max",
    "final_rel_residual_max","sec_per_work_unit_from_campaign",
    "estimated_cumulative_pressure_time_s"
]
with (OUT/"pressure_history_window_summary.csv").open("w",newline="") as fh:
    wr=csv.DictWriter(fh,fieldnames=fields); wr.writeheader(); wr.writerows(all_rows)

# Compact solver comparison: if multiple files map to same solver, take the first
# one with the most rows (expected one file per solver for the requested condition).
best={}
for solver,p,n,cc in selected_files:
    if solver=="UNKNOWN": continue
    if solver not in best or n>best[solver][1]:
        best[solver]=(p,n)
compact=[]
for solver,(p,n) in best.items():
    subset=[r for r in all_rows if r["solver"]==solver and r["source_file"]==str(p)]
    compact.extend(subset)
with (OUT/"pressure_history_compact.csv").open("w",newline="") as fh:
    wr=csv.DictWriter(fh,fieldnames=fields); wr.writeheader(); wr.writerows(sorted(compact,key=lambda r:(r["window_end_s"],r["solver"])))

# Crossover estimates based on average campaign cost/work-unit.
# Compare cumulative estimated pressure time at each common window.
bywin={}
for r in compact:
    bywin.setdefault(r["window_end_s"],{})[r["solver"]]=r
with (OUT/"pressure_history_estimated_crossover.txt").open("w") as fh:
    fh.write("Estimated cumulative pressure-time comparison\n")
    fh.write("Uses campaign-average seconds/work-unit, not per-step measured timing.\n\n")
    for w in WINDOWS:
        d=bywin.get(w,{})
        if not d: continue
        fh.write(f"t <= {w:g} s\n")
        vals=[]
        for sol,r in d.items():
            vals.append((r["estimated_cumulative_pressure_time_s"],sol,r))
        for est,sol,r in sorted(vals):
            fh.write(f"  {sol:5s}: work={r['cumulative_work_units']:.6g}, mean={r['mean_work_units_per_solve']:.6g}, zero={100*r['zero_work_fraction']:.2f}%, estPressure={est:.6g} s\n")
        fh.write("\n")

print("ANALYSIS_COMPLETE")
print("files",len(files),"selected",len(selected_files),"solvers",sorted(best))
