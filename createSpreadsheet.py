#!/usr/bin/env python

## Creates spreadsheets combining input metadata and noteworthy results
import sys
import pandas as pd
from pathlib import Path
import json

def getMetaDF(meta_file):
    """Returns metadata as df"""

    df = pd.read_csv(meta_file)
    # finalize some column names
    newCols = {"Seq ID":"Sample #","Ct N1":"Ct N gene","Qubit reading":"post-ARTIC Qubit (ng/ul)","Zip Code":"Zip code","Sample date":"Test date"}
    df = df.rename(columns=newCols)
    cols = df.columns
    if "Zip code" in cols:
        df["Zip code"] = df["Zip code"].astype('Int64')
    # make dates dates (to ensure they're consistent)
    for col in cols:
        if "date" in col.lower():
            df[col] = pd.to_datetime(df[col])
    print(f"metadata\t{len(df)}")
    return df

def get_ambig(x):
    if x == 'failed to map': return 1.0
    try:
        return float(x.split(":")[-1])
    except(ValueError):
        print(f"'Ambiguity_content' column contains unaccounted for string '{x}'")
        exit(1)

def getPangoDF(pango_file):
    """Returns useful pangolin metadata as df"""

    df = pd.read_csv(pango_file)
    df["status"] = df["qc_status"].apply(lambda x: x.replace("pass","passed_qc"))
    pangolin_cols = ["taxon","lineage","scorpio_call","ambiguity_score","is_designated","status","qc_notes","note","pangolin_version"]
    df = df[pangolin_cols]
    df["Sample #"] = df.taxon.apply(lambda x: x.split("/")[0])
    df = df.rename(columns={
        "scorpio_call":"Scorpio lineage",
        "lineage":"Pango lineage",
        "status":"Pangolin QC",
        "ambiguity_score":"Scorp Score",
        "note":"pango_note",
        "qc_notes":"Ambiguous_content"})
    df["Ambiguous_content"] = df["Ambiguous_content"].apply(lambda x: get_ambig(x)) # convert from "Ambiguous_content:0.02" --> 0.02
    print(f"pangolin\t{len(df)}")
    print(df)
    return df

def getNextcladeDF(nextclade_file,nextclade_version):
    """Returns nextclade file as df with useful columns"""

    df = pd.read_csv(nextclade_file,sep=";")
    df = df[["seqName","clade","qc.overallScore","qc.overallStatus","totalMissing","aaSubstitutions","substitutions"]]
    df["Sample #"] = df["seqName"].apply(lambda x: x.split("/")[0])
    df = df.rename(columns={
        "clade":"Nextstrain Clade",
        "qc.overallScore":"Nextstrain QC",
        "qc.overallStatus":"Nextstrain Status",
        "totalMissing":"Total Missing",
        "aaSubstitutions":"AA Substitutions",
        "substitutions":"Nucleotide Substitutions"})
    df = df.drop_duplicates()
    if "Total Missing" in df.columns:
        df["Total Missing"] = df["Total Missing"].astype('Int64')
    df["nextclade_version"] = nextclade_version
    print(f"nextclade\t{len(df)}")
    print(df)
    return df

def getReadsCountsDF(reads_file,counts="Reads on Barcode"):
    """Returns read counts by barcode as df"""

    df = pd.read_csv(reads_file,header=None,names=["Sample #",counts])
    print(f"reads\t\t{len(df)}")
    print(df)
    return df

def checkDfs(dfs):
    """Verifies dfs aren't empty; warns if they're awfully small"""

    for name,df in dfs.items():
        if df.empty:
            raise Exception(f"No data found in {name} df.")

def writeSpreadsheet(df,outdir,plate,which_report="All"):
    """Writes a spreadsheet for individual source labs"""

    # write out df if not empty
    if len(df) > 0:
        print("outdir:",outdir)
        outfile = outdir / f"Sequencing-report-{plate}-{which_report}.csv"
        print(f"Writing out {which_report} report\n\tOutfile: {outfile}")
        df.to_csv(outfile,index=False)
    else: print(f"{which_report}\n\tSkipping. No samples found in this batch for Source Lab: '{which_report}'")
    
def writeSequencingReport(combo,outdir,plate,reportMap):
    """Writes df to csv"""

    # output file to csvs (split by source) and edit permissions
    print(f"\nWriting out spreadsheet CSVs:")
    # this stuff vvv is mirrored in Prepare-Corvaseq.py -- ensure it stays updated together
    # reportMap = {"Campus":["Campus",None],"Mecklenburg":["Starmed|StarMed","MCPH"],"StarMed":["Starmed|StarMed",None],"Mission":["Mission",None]}            # EDIT if more source labs added

    # if there are multiple source labs, write a report file for each
    if "Source Lab" in combo.columns and len(combo["Source Lab"].unique()) > 1:
        if reportMap == '':
            # skip individual reports
            pass
        elif reportMap == "{}":
            # create reports for each source lab
            for source_lab in combo["Source Lab"].unique():
                df = combo[combo["Source Lab"].str.contains(source_lab, na=False)]
                writeSpreadsheet(df,outdir,plate)
        else:
            # create reports for each group specified in reportMap
            reportMap = json.loads(reportMap.replace("'",'"'))
            for which_report, data in reportMap.items():
                source_lab = data.get("source_lab",which_report)
                limit_to = data.get("limit_to",None)
                df = combo[combo["Source Lab"].str.contains(source_lab, na=False)]
                if limit_to:
                    df = combo[combo["Sample #"].str.contains(limit_to, na=False)]
                writeSpreadsheet(df,outdir,plate,which_report)

    # write out cumulative report (All)
    outfile = outdir / f"Sequencing-report-{plate}-All.csv"
    print(f"Writing out combined report\n\tOutfile: {outfile}")
    combo.to_csv(outfile, index=False)

def main():
    meta_file = Path(sys.argv[1])
    pangolin_file = Path(sys.argv[2])
    nextclade_file = Path(sys.argv[3])
    nextclade_version = Path(sys.argv[4])
    read_counts_file = Path(sys.argv[5])
    raw_read_counts_file = Path(sys.argv[6])
    plate = sys.argv[7]
    outdir = Path(sys.argv[8])
    reportMap = sys.argv[9]

    print("Inputs:")
    print("meta_file:",meta_file)
    print("pangolin_file:",pangolin_file)
    print("nextclade_file:",nextclade_file)
    print("nextclade_version:",nextclade_version)
    print("read_counts_file:",read_counts_file)
    print("raw_read_counts_file:",raw_read_counts_file)
    print("plate:",plate)
    print("outdir:",outdir)
    print("reportMap:",reportMap)
    print()

    # read in all csvs to dataframes
    print("Reading in CSVs\ndf\t\tlength")
    meta_df = getMetaDF(meta_file)
    pangolin_df = getPangoDF(pangolin_file)
    nextclade_df = getNextcladeDF(nextclade_file,nextclade_version)
    reads_df = getReadsCountsDF(read_counts_file)
    raw_reads_df = getReadsCountsDF(raw_read_counts_file,"read_counts_raw")

    # ensure all dfs have data
    dfs = {"Metadata":meta_df,"Nextclade":nextclade_df,"Pangolin":pangolin_df,"read counts":reads_df,"raw read counts":raw_reads_df}
    checkDfs(dfs)
    # for df in dfs.values():
    #     df["Sample #"] = df["Sample #"].astype(str)

    # merge dfs
    combo:pd.DataFrame = meta_df.merge(pangolin_df,on="Sample #",how="outer"
                ).merge(nextclade_df,on="Sample #",how="outer"
                ).merge(reads_df,on="Sample #",how="outer"
                ).merge(raw_reads_df,on="Sample #",how="outer")

    # write out results
    writeSequencingReport(combo,outdir,plate,reportMap)
    # combo.to_csv(outfile,index=False)

    print(f"\ndone\n")

if __name__ == "__main__":
    main()
