#!/usr/bin/env python
# -*- coding: utf-8 -*-

import multiprocessing
import numpy as np
import pandas as pd
from Library.spectrum_library import Database
from simRank_library import *
import os
import time
import networkx as nx
from typing import *


def simRank_network(
    sample_folder,
    output_dir,
    scorematrix_path,
    ms2delta,
    ppm,
    maxrank,
    simrank_thr,
    peaktable_path=None,
    n_cpu=None,
    pm_tolerance=1e15,
    msdelta=0.01,
):
    """
    filter one sample file based on simRank
    :param sample_folder: folder containing sample file
    :param peaktable_path: location of peaktable
    :param output_dir: output directory
    :param scorematrix_path: location of scorematrix for simRank calculation
    :return: graphml file
    """
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

    # generate rankmz dataframe of sample file
    count = 0
    for file_name in os.listdir(sample_folder):
        if file_name.lower().endswith(".mzxml"):
            count += 1
    if count < 1:
        print("at least one mzxml file should be in the folder")
        exit()

    db = Database(sample_folder, os.path.basename(sample_folder))
    db.load_files(
        inten_thresh=inten_thresh,
        engine=engine,
        n_workers=n_workers,
        remove_pre=remove_precursor,
        min_mz_num=min_mz_num,
    )
    db.merge_spectra(
        ms2delta=ms2delta, ppm=ppm, rt=rt, by_col_e=if_merge_samples_byenergy
    )
    db.normalize_peaks()
    db_frame = get_database_rankmz(db, maxrank=maxrank)
    db_frame = db_frame.reset_index(drop=True)
    # filter entries in db_frame based on peaktable
    if (peaktable is not None) and (peaktable is not False):
        if len(db_frame) > 0:
            index = []
            mzlist = np.array(peaktable["mz"])
            rtlist = np.array(peaktable["rt"])
            for i in range(len(db_frame)):
                mz_i = db_frame.loc[i, "precursor"]
                rt_i = db_frame.loc[i, "rt"]
                tol = abs(round(min(mz_i * ppm * 1e-6, msdelta), 4))
                matchlist = (abs(mzlist - mz_i) < tol) & (abs(rtlist - rt_i) < rt)
                if True in matchlist:
                    index.append(i)
                    try:
                        temp = (peaktable.loc[matchlist]).reset_index()
                        db_frame.loc[i, "inten"] = temp.loc[
                            0, "inten"
                        ]  # if peaktable contains intensity info, use that value
                    except:
                        print(
                            "no intensity given in peaktable, intensity will be based on MS2 precursor intensity"
                        )
            db_frame = db_frame.loc[index]
            db_frame = db_frame.reset_index(drop=True)

    # graphml file cannot take NoneType values
    db_frame = db_frame.fillna(value="NA")
    # generate network
    nodeinfo = db_frame[
        ["ID", "precursor", "inten", "rt", "iontype", "charge"]
    ].to_dict("index")
    a = time.time()
    G = nx.Graph()
    G.add_nodes_from(nodeinfo)
    nx.set_node_attributes(G, nodeinfo)

    edge_list = returnEdgeMult(
        db_frame, ms2delta, ppm, maxrank, scorematrix, simrank_thr, n_cpu, pm_tolerance
    )
    G.add_weighted_edges_from(edge_list)
    print("time for generating network: ", time.time() - a)

    # write graphml file
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    nx.write_graphml(G, os.path.join(output_dir, "simrank_network_mult.graphml"))

    print(time.time())
    print(time.time() - a)


def returnEdgeWithWeightFromTwoFrames(
    db_frame1, db_frame2, ms2delta, ppm, maxrank, scorematrix, simrank_thr, pm_tolerance
):
    result = []
    for i in db_frame1.index:
        for j in db_frame2.index:
            if i == j:
                continue
            rankmz1 = db_frame1.loc[i, "rankmz"]
            rankmz2 = db_frame2.loc[j, "rankmz"]
            precursor1 = db_frame1.loc[i, "precursor"]
            precursor2 = db_frame2.loc[j, "precursor"]
            if abs(precursor1 - precursor2) > pm_tolerance:
                continue
            simrankmat1 = get_simrankmat_pre_mod(
                rankmz1, rankmz2, precursor1, precursor2, ms2delta, ppm, maxrank
            )
            simrankmat2 = get_simrankmat_pre_mod(
                rankmz2, rankmz1, precursor2, precursor1, ms2delta, ppm, maxrank
            )
            score = max(
                round(np.sum(simrankmat1 * scorematrix), 4),
                round(np.sum(simrankmat2 * scorematrix), 4),
            )
            if (
                score >= simrank_thr
            ):  # if two nodes simrank similar, then generate an edge
                result.append((i, j, score))
    return result


def returnEdgeWithWeight(
    db_frame, i, j, ms2delta, ppm, maxrank, scorematrix, simrank_thr, pm_tolerance
):
    result: List[set] = [()]
    try:
        db_frame_2 = db_frame.loc[i : j - 1]
        result = returnEdgeWithWeightFromTwoFrames(
            db_frame1=db_frame,
            db_frame2=db_frame_2,
            ms2delta=ms2delta,
            ppm=ppm,
            maxrank=maxrank,
            scorematrix=scorematrix,
            simrank_thr=simrank_thr,
            pm_tolerance=pm_tolerance,
        )
    except Exception as e:
        print(e)
    return result


def singleArg_returnEdgeWithWeight(Arg: dict):
    db_frame = Arg["db_frame"]
    i = Arg["start"]
    j = Arg["end"]
    ms2delta = Arg["ms2delta"]
    ppm = Arg["ppm"]
    maxrank = Arg["maxrank"]
    scorematrix = Arg["scorematrix"]
    simrank_thr = Arg["simrank_thr"]
    pm_tolerance = Arg["pm_tolerance"]
    result: List[set] = returnEdgeWithWeight(
        db_frame, i, j, ms2delta, ppm, maxrank, scorematrix, simrank_thr, pm_tolerance
    )
    return result


def returnEdgeMult(
    db_frame, ms2delta, ppm, maxrank, scorematrix, simrank_thr, n_cpu, pm_tolerance
):
    """
    multprocess function for generating edges
    """
    if n_cpu is None:
        n_cpu = int(multiprocessing.cpu_count())
    pool = multiprocessing.Pool(processes=n_cpu)
    db_frame_length = len(db_frame)

    if db_frame_length > 1:
        split_interval = int(db_frame_length / n_cpu)
        if split_interval == 0:
            split_index = [0, db_frame_length]
        else:
            split_index = [i * split_interval for i in range(n_cpu)]
            split_index.append(db_frame_length)
        split_start = split_index[:-1]
        split_end = split_index[1:]
        edge_list = []

        args_list: list = []
        for i in range(len(split_start)):
            start = split_start[i]
            end = split_end[i]
            args_list.append(
                {
                    "db_frame": db_frame,
                    "start": start,
                    "end": end,
                    "ms2delta": ms2delta,
                    "ppm": ppm,
                    "maxrank": maxrank,
                    "scorematrix": scorematrix,
                    "simrank_thr": simrank_thr,
                    "pm_tolerance": pm_tolerance,
                }
            )

        for i in range(len(split_start)):
            edge_list.append(
                pool.apply_async(
                    singleArg_returnEdgeWithWeight,
                    (args_list[i],),
                )
            )
        pool.close()
        pool.join()
        edge_list = [item.get() for item in edge_list]
        edge_list = [item for sublist in edge_list for item in sublist]
        edge_list = [item for item in edge_list if len(item) != 0]
    else:
        print("not enough nodes to generate network")

    return edge_list


if __name__ == "__main__":
    sample_folder = r""
    peaktable_path = r""
    output_dir = r""
    scorematrix_path = r""
    inten_thresh = 1
    engine = r"utf-8"
    n_workers = multiprocessing.cpu_count()
    rt = 30
    ppm = 50
    pm_tolerance = 300
    msdelta = 0.01
    ms2delta = 0.01
    min_mz_num = 1
    if_merge_database = True
    if_merge_database_byenergy = True
    if_merge_samples_byenergy = False
    remove_precursor = True
    hit_num = 5
    maxrank = 30
    simrank_thr = 16

    simRank_network(
        sample_folder=sample_folder,
        output_dir=output_dir,
        scorematrix_path=scorematrix_path,
        ms2delta=ms2delta,
        ppm=ppm,
        maxrank=maxrank,
        simrank_thr=simrank_thr,
        peaktable_path=peaktable_path,
        n_cpu=n_workers,
        pm_tolerance=pm_tolerance,
        msdelta=msdelta,
    )
    print("done")
