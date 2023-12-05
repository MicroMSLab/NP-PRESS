#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
import numpy as np
import pandas as pd
from Library.spectrum_library import Database, SpectrumCollection
from simRank_library import *
import os
import datetime


def simRank_filter(
    sample_path, control_df, output_dir, scorematrix_path, peaktable_path=None
):
    """
    filter one sample file based on simRank
    :param sample_path: location of sample file
    :param peaktable_path: location of peaktable
    :param control_df: delaed reference/control folder
    :param output_dir: output directory
    :param scorematrix_path: location of scorematrix for simRank calculation
    :return: no return value, output simRank results and filtered peaktable to output directory
    """

    # check output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # load scorematrix
    if os.path.exists(scorematrix_path) and scorematrix_path.lower().endswith(".csv"):
        scorematrix = np.array(pd.read_csv(scorematrix_path))
    else:
        print("scorematrix is not given")
        return

    # load peaktable of peaks for analysis
    if os.path.exists(peaktable_path) and peaktable_path.lower().endswith(".csv"):
        peaktable = pd.read_csv(peaktable_path)
        peaktable = peaktable.reset_index(drop=True)
        if (
            ("mz" in peaktable.columns)
            and ("rt" in peaktable.columns)
            and (len(peaktable) > 0)
        ):
            try:
                peaktable = peaktable.loc[
                    :, ["mz", "rt", "inten"]
                ]  # if peaktable contains intensity info, include it
            except:
                peaktable = peaktable.loc[:, ["mz", "rt"]]
            print("Peaktable is loaded!")
        else:
            print("Peaktable is not given!")
            peaktable = False
    else:
        print("Peaktable is not given!")
        peaktable = False

    # load sample file, must be .mzxml file
    if os.path.exists(sample_path) and sample_path.lower().endswith(".mzxml"):
        query_collection = SpectrumCollection(sample_path, type="ms2")
        query_collection.load_from_file(
            inten_thresh=inten_thresh,
            peaktable=peaktable,
            drop_mslevels=[1],
            engine=engine,
            remove_pre=remove_precursor,
            min_mz_num=min_mz_num,
        )
        query_collection.make_ms1list_sequential()
        query_collection.merge_spectra(
            ppm=ppm, rt=rt, ms2delta=ms2delta, by_col_e=if_merge_samples_byenergy
        )
        query_collection.normalize_peaks()
    else:
        print("sample file is not correct")
        return

    # load database Dataframe files, match against each database dataframe

    result = simrank_speccol2db(
        query_collection,
        control_df,
        threshold=simrank_thr,
        scorematrix=scorematrix,
        hit_num=hit_num,
        pre_diff=pm_tolerance,
        maxrank=maxrank,
        delta=ms2delta,
        ppm=ppm,
        n_cpu=n_workers,
    )

    result = result.reset_index(drop=True)

    # output new peaktable after filtering
    if (peaktable is not None) and (peaktable is not False):
        if len(result) > 0:
            index = []
            mzlist = np.array(result["signal_mz"])
            rtlist = np.array(result["signal_rt"])
            for i in range(len(peaktable)):
                mz_i = peaktable.loc[i, "mz"]
                rt_i = peaktable.loc[i, "rt"]
                tol = abs(round(min(mz_i * ppm * 1e-6, msdelta), 4))
                matchlist = (abs(mzlist - mz_i) < tol) & (abs(rtlist - rt_i) < rt)
                if True not in matchlist:
                    index.append(i)
            peaktable = peaktable.loc[index, :]
        output_file = os.path.join(
            output_dir,
            "simRank_filteredPeaktable_"
            + str(datetime.datetime.now().strftime("%Y%m%d%H%M%S"))
            + ".csv",
        )
        peaktable.to_csv(output_file, index=False)


def Deal_Control(control_folder_dir_list, if_merge=True):
    """
    load database files, if more than one database folder than combine them, generate pickle files
    :param db_dir_list:
    :param if_merge:
    :return:
    """
    db = Database("", "")
    for folder in control_folder_dir_list:
        foldername = os.path.basename(folder)
        if foldername.startswith("."):
            continue
        db_new = Database(folder, foldername)
        db_new.load_files(
            inten_thresh=inten_thresh,
            engine=engine,
            n_workers=n_workers,
            remove_pre=remove_precursor,
            min_mz_num=min_mz_num,
        )
        db_new.normalize_peaks()
        db = db.combine_database(db_new)
        del db_new

    if if_merge:
        db.merge_spectra(ms2delta=ms2delta)
        db.normalize_peaks()

    db_frame = get_database_rankmz(db, maxrank=maxrank)
    del db

    return db_frame


if __name__ == "__main__":
    # set sample and database location
    sample_path = r".mzXML"
    peaktable_path = (
        r""  # you can include 'inten' in your peaktable, if you only have one sample
    )
    control_folder_dir_list = [r"/control"]
    output_dir = "/result"
    scorematrix_path = "../scorematrix.csv"
    inten_thresh = 1
    engine = "utf-8"
    n_workers = multiprocessing.cpu_count()
    rt = 30
    ppm = 20
    pm_tolerance = 300
    msdelta = 0.01
    ms2delta = 0.01
    min_mz_num = 1
    if_merge_database = True
    if_merge_database_byenergy = (
        True  # this parameter only works for database generated from .mzXML
    )
    if_merge_samples_byenergy = False
    remove_precursor = True
    hit_num = 5
    maxrank = 30
    simrank_thr = 0

    simRank_filter(
        sample_path=sample_path,
        control_df=Deal_Control(
            control_folder_dir_list=control_folder_dir_list, if_merge=if_merge_database
        ),
        output_dir=output_dir,
        scorematrix_path=scorematrix_path,
        peaktable_path=peaktable_path,
    )
    print("done")
