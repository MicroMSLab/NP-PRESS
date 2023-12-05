#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from Library.spectrum_library import Database, SpectrumCollection
import copy, os, scipy, multiprocessing, math


# TODO
# 对于一个database进行rankmz计算 done
# 对于一个spec计算它和database中所有对象的simrank_score done
# 对于一个database计算每个pair的rankmz的匹配情况
# 根据rankmz的匹配情况和分子指纹的匹配程度计算H0matrix和H1matrix

def gammapdf_fit(x, y, p0=None, maxfev=500000):
    def gamma_pdf(x, A, B, C):
        return B ** A / (scipy.special.gamma(A)) * x ** (A - 1) * math.e ** (-B * x) * C

    A, B, C = scipy.optimize.curve_fit(gamma_pdf, x, y, p0=p0, maxfev=maxfev)[0]
    ypred = B ** A / (scipy.special.gamma(A)) * x ** (A - 1) * math.e ** (-B * x) * C
    return A, B, C, ypred

def linear_fit(x, y, maxfev=500000):
    def f(x, A, B):
        return A * x + B
    A, B = scipy.optimize.curve_fit(f, x, y, maxfev=maxfev)[0]
    ypred = A * x + B
    return ypred

def calTaniScore(FP1, FP2):
    try:
        return(float((str(bin(int(FP1)&int(FP2))).count('1'))/(str(bin(int(FP1)|int(FP2))).count('1'))))
    except ValueError:
        # print('{} or {} is not decimal fingerprints'.format(str(FP1),str(FP2)))
        return -1

def calCosineScore(FP1, FP2):
    import math
    try:
        return (float((str(bin(int(FP1) & int(FP2))).count('1')) /
                      (math.sqrt(str(bin(int(FP1) | int(0))).count('1'))) /
                      (math.sqrt(str(bin(int(FP2) | int(0))).count('1')))
                      )
                )
    except ValueError:
        # print('{} or {} is not decimal fingerprints'.format(str(FP1),str(FP2)))
        return -1

def get_rankmz(peaks):
    peaks = np.asarray(peaks)
    order = np.argsort(peaks[:, 1])
    order = order[::-1]
    return peaks[order, 0]


def get_database_rankmz(db, maxrank=30):
    '''
    get rankmz array for database, the database can be merged or not merged
    :param db: Database object
    :param maxrank: maximum rank, e.g. when maxrank=30, only the top 30 peaks will be kept
    :return: pandas DataFrame with columns=['ID', 'precursor', 'inten', 'rt', 'iontype', 'charge',  'rankmz']
    '''
    resultlist = [pd.DataFrame()]
    for speccol in db.SpectrumCollection_list:
        ID = speccol.MID
        for spectrum in speccol.spectrum_list:
            precursor = spectrum.mz
            iontype = spectrum.ion
            charge = spectrum.charge
            collision_energy = spectrum.collision_energy
            rankmz = get_rankmz(spectrum.peaks)
            rt = spectrum.retention_time
            inten = spectrum.precursor_intensity #use precursor intensity as the signal intensity
            if len(rankmz) > maxrank:
                rankmz = rankmz[0:maxrank]
            if (speccol.MID is None) or (speccol.MID == 'None'):#if generated from mzxml, then no MID is given
                ID = os.path.basename(spectrum.filename).split('.')[0]+\
                     '_m/z:'+str(round(float(spectrum.mz),2))\
                     +'_'+str(round(float(spectrum.retention_time),1))+'s'
            # resultlist.append(
            #     pd.DataFrame([[ID, precursor, iontype, charge, collision_energy, rankmz, spectrum]],
            #                  columns=['ID', 'precursor', 'iontype', 'charge', 'collision_energy', 'rankmz', 'spec']))
            resultlist.append(
                pd.DataFrame([[ID, round(precursor,4), round(inten, 2), round(rt,2), iontype, charge, collision_energy, rankmz]],
                             columns=['ID', 'precursor', 'inten', 'rt', 'iontype', 'charge', 'collision_energy', 'rankmz']))
    return pd.concat(resultlist, ignore_index=True)


def get_simrankmat_pre(rankmz1, rankmz2, precursor1, precursor2, delta=0.01, ppm=20, maxrank=30):
    '''
    given rankmz of two compounds and their precursor mz, calculate the rank similarity matrix
    :param rankmz1: a list of ordered mz according to their intensity from high to low
    :param rankmz2: a list of ordered mz according to their intensity from high to low
    :param precursor1: precursor mz of the first compound
    :param precursor2: precursor mz of the second compound
    :param delta: match tolerance for MS2 peaks
    :param ppm: match tolerance for MS2 peaks
    :param maxrank: maximum rank, e.g. when maxrank=30, only the top 30 peaks will be kept
    :return: a 0-1 matrix with dimension of (maxrank+1)*maxrank
    '''
    if len(rankmz1) > maxrank:
        rankmz1 = rankmz1[0:maxrank]
    if len(rankmz2) > maxrank:
        rankmz2 = rankmz2[0:maxrank]

    result = np.zeros((maxrank + 1, maxrank))
    # if rank1 of rankmz1 matches rank2 of rankmz2, then result[rank1, rank2] = 1
    # if rank1 of rankmz1 does not match any of rankmz2, then result[maxrank, rank1] = 1
    for rank1 in range(len(rankmz1)):
        mz1 = rankmz1[rank1]
        tolerance = min(delta, mz1 * ppm * 1e-6)  # tolerance is determined based on delta and ppm
        match = False
        for rank2 in range(len(rankmz2)):
            mz2 = rankmz2[rank2]
            if (abs(mz1 - mz2) < tolerance) or (abs(mz1 - mz2 - (precursor1 - precursor2)) < tolerance):
                result[rank1, rank2] = 1
                match = True
                break
        if not match:
            result[maxrank, rank1] = 1
    return result

def get_simrankmat_pre_mod(rankmz1, rankmz2, precursor1, precursor2, delta=0.01, ppm=20, maxrank=30):
    '''
    given rankmz of two compounds and their precursor mz, calculate the rank similarity matrix, penalty based on matched mz, lower mz will have penalty
    :param rankmz1: a list of ordered mz according to their intensity from high to low
    :param rankmz2: a list of ordered mz according to their intensity from high to low
    :param precursor1: precursor mz of the first compound
    :param precursor2: precursor mz of the second compound
    :param delta: match tolerance for MS2 peaks
    :param ppm: match tolerance for MS2 peaks
    :param maxrank: maximum rank, e.g. when maxrank=30, only the top 30 peaks will be kept
    :return: a 0-1 matrix with dimension of (maxrank+1)*maxrank
    '''
    if len(rankmz1) > maxrank:
        rankmz1 = rankmz1[0:maxrank]
    if len(rankmz2) > maxrank:
        rankmz2 = rankmz2[0:maxrank]

    result = np.zeros((maxrank + 1, maxrank))
    # if rank1 of rankmz1 matches rank2 of rankmz2, then result[rank1, rank2] = 1
    # if rank1 of rankmz1 does not match any of rankmz2, then result[maxrank, rank1] = 1
    for rank1 in range(len(rankmz1)):
        mz1 = rankmz1[rank1]
        tolerance = min(delta, mz1 * ppm * 1e-6)  # tolerance is determined based on delta and ppm
        match = False
        for rank2 in range(len(rankmz2)):
            mz2 = rankmz2[rank2]
            if (abs(mz1 - mz2) < tolerance) or (abs(mz1 - mz2 - (precursor1 - precursor2)) < tolerance):
                # mz match lower than 100 will have penalty
                result[rank1, rank2] = min(1, round(min(mz1,mz2)/100, 2))
                match = True
                break
        if not match:
            result[maxrank, rank1] = 1
    return result

def get_simrankmat(rankmz1, rankmz2, precursor1=0, precursor2=0, delta=0.01, ppm=20, maxrank=30):
    '''
    given rankmz of two compounds, calculate the rank similarity matrix, not considering precursor difference
    :param rankmz1: a list of ordered mz according to their intensity from high to low
    :param rankmz2: a list of ordered mz according to their intensity from high to low
    :param precursor1: precursor mz of the first compound
    :param precursor2: precursor mz of the second compound
    :param delta: match tolerance for MS2 peaks
    :param ppm: match tolerance for MS2 peaks
    :param maxrank: maximum rank, e.g. when maxrank=30, only the top 30 peaks will be kept
    :return: a 0-1 matrix with dimension of (maxrank+1)*maxrank
    '''
    if len(rankmz1) > maxrank:
        rankmz1 = rankmz1[0:maxrank]
    if len(rankmz2) > maxrank:
        rankmz2 = rankmz2[0:maxrank]

    result = np.zeros((maxrank + 1, maxrank))
    # if rank1 of rankmz1 matches rank2 of rankmz2, then result[rank1, rank2] = 1
    # if rank1 of rankmz1 does not match any of rankmz2, then result[maxrank, rank1] = 1
    for rank1 in range(len(rankmz1)):
        mz1 = rankmz1[rank1]
        tolerance = min(delta, mz1 * ppm * 1e-6)  # tolerance is determined based on delta and ppm
        match = False
        for rank2 in range(len(rankmz2)):
            mz2 = rankmz2[rank2]
            if (abs(mz1 - mz2) < tolerance):
                result[rank1, rank2] = 1
                match = True
                break
        if not match:
            result[maxrank, rank1] = 1
    return result

#TODO: use deeplearning model to calculate score
def get_simrankscore(rankmz1, precursor1, db_frame, scorematrix, pre_diff=300, maxrank=30, delta=0.01, ppm=20, cpuid=None):
    '''
    called in get_ranksimhit_multi(), calculate ranksim hits in mono-processor mode
    :param rankmz1:
    :param precursor1:
    :param db_frame:
    :param scorematrix:
    :param delta:
    :param ppm:
    :param cpuid:
    :return:
    '''
    db_frame = copy.deepcopy(db_frame.reset_index(drop=True))
    for i in range(len(db_frame)):
        precursor2 = db_frame.loc[i, 'precursor']
        rankmz2 = db_frame.loc[i, 'rankmz']
        # difference between the precursors of the two compounds must be smaller than threshold
        if abs(precursor1 - precursor2) > abs(pre_diff):
            continue
        # calculate the ranksim matrix between the spec and the database compound
        simrankmat1 = get_simrankmat_pre(rankmz1, rankmz2, precursor1, precursor2, delta, ppm, maxrank)
        simrankmat2 = get_simrankmat_pre(rankmz2, rankmz1, precursor2, precursor1, delta, ppm, maxrank)
        db_frame.loc[i, 'simrank_score'] = max(round(np.sum(simrankmat1 * scorematrix), 4),round(np.sum(simrankmat2 * scorematrix), 4))
        if (i+1) % 1000 == 0:
            print('\rprocessor {}, {} compounds compared'.format(cpuid, i+1), end='', flush=True)
    columns = list(db_frame.columns)
    columns.remove('rankmz')
    db_frame = db_frame.loc[:, columns]
    return db_frame

#TODO: use deeplearning model to calculate score
def get_simrankscore_multi(spec, db_frame, scorematrix, hit_num, pre_diff=300, maxrank=30, delta=0.01, ppm=20,
                           n_cpu=None):
    '''
    compare a spectrum against a rankmz-calculated database using ranksim algorithm, based on the given
    scorematrix
    :param spec: spectrum to search against the database
    :param db_frame: db DataFrame containing compound info and rankmz
    :param scorematrix: scorematrix [(maxrank+1)*maxrank] to calculate similarity
    :param hit_num: how many hits to report
    :param pre_diff: precursor difference tolerance
    :param maxrank: maximum rank of peaks to consider
    :param delta: MS2 peak tolerance
    :param ppm: MS2 peak tolerance
    :param n_cpu: how many cpus to use in parellel
    :return: DataFrame with hit info and similarity score, columns=['ID', 'precursor', 'iontype', 'charge',  'simrank_score']
    '''
    if len(db_frame) < 1:
        print('warning: database is empty!')
        return
    precursor1 = spec.mz
    rankmz1 = get_rankmz(spec.peaks)
    if len(rankmz1) > maxrank:
        rankmz1 = rankmz1[0:maxrank]
    db_frame = db_frame.reset_index(drop=True)
    # multiprocessing
    totalnum = len(db_frame)
    pool = multiprocessing.Pool(processes=n_cpu)
    if n_cpu is None:
        n_cpu = multiprocessing.cpu_count()
    Ms = []
    base = int(totalnum / n_cpu)
    if totalnum < n_cpu:
        for i in range(totalnum):
            Ms.append(pool.apply_async(
                get_simrankscore, (rankmz1, precursor1,
                                   db_frame.loc[i, :],
                                   scorematrix, pre_diff, maxrank, delta, ppm, i + 1)
            ))
    else:
        for i in range(n_cpu):
            if i+1 < n_cpu:
                # .loc include the start and end index
                Ms.append(pool.apply_async(
                    get_simrankscore,(rankmz1, precursor1,
                                     db_frame.loc[i * base:
                                                           (i + 1) * base - 1, :],
                                     scorematrix, pre_diff, maxrank, delta, ppm, i + 1)
                ))
            else:
                # .loc include the start and end index
                Ms.append(pool.apply_async(
                    get_simrankscore,(rankmz1, precursor1,
                                     db_frame.loc[i * base:
                                                           totalnum, :],
                                     scorematrix, pre_diff, maxrank, delta, ppm, i + 1)
                ))
    pool.close()
    pool.join()
    Ms = [M.get() for M in Ms]
    hit_table = pd.concat(Ms, ignore_index=True)
    # sort by simrank_score, return the top hits
    hit_table = hit_table.sort_values(by='simrank_score', ascending=False)
    hit_table = hit_table.reset_index(drop=True)
    if len(hit_table) > hit_num:
        hit_table = hit_table.loc[0:hit_num, :]
    return hit_table

#TODO: use deeplearning model to calculate score
def simrank_speccol2db(speccol, db_frame, threshold, scorematrix, hit_num, pre_diff=300, maxrank=30, delta=0.01, ppm=20,
                       n_cpu=None):
    result = [pd.DataFrame()]
    for spec in speccol.spectrum_list:
        temp = get_simrankscore_multi(spec, db_frame, scorematrix=scorematrix, hit_num=hit_num, pre_diff=pre_diff,
                                      maxrank=maxrank, delta=delta, ppm=ppm, n_cpu=n_cpu)
        temp['file'] = os.path.basename(spec.filename).split('.')[0]
        temp['signal_mz'] = spec.mz
        temp['signal_rt'] = spec.retention_time
        temp = temp.loc[temp['simrank_score'] >= threshold, :]
        if len(temp) > 0:
            result.append(temp)
    result = pd.concat(result, ignore_index=True)
    result = result.reset_index(drop=True)
    return result


if __name__ == '__main__':
    # set sample and database location
    sample_path = ''
    control_path = ''
    db_dir_list = ['']
    output_path = ''
    scorematrix_path = ''
    inten_thresh = 1
    engine = 'utf-8'
    n_workers = multiprocessing.cpu_count()
    rt = 30
    ppm = 20
    pm_tolerance = 300
    ms2delta = 0.01
    min_matched_peaks = 3
    min_score = 0.0
    MID = None
    min_mz_num = 1
    if_merge_database = True
    if_merge_database_byenergy = True  # this parameter only works for database generated from .mzXML
    if_merge_samples_byenergy = True
    remove_precursor = True
    hit_num = 5
    maxrank = 30
    simrank_thr = 0
    
    # get database frame containing ranked mz
    for db_dir in db_dir_list:
        if db_dir.startswith('.'):
            continue
        if not (os.path.isdir(db_dir)):
            print('warning: database folder {} is not correct'.format(db_dir))
            continue
        try:
            db = Database(db_dir, os.path.basename(db_dir))
            db.load_files(inten_thresh=inten_thresh, engine=engine,
                          n_workers=n_workers, remove_pre=remove_precursor, min_mz_num=min_mz_num)
        except:
            print('warning: database folder {} cannot be loaded'.format(db_dir))
            continue
        if if_merge_database:
            db.merge_spectra(ppm=ppm, rt=rt, ms2delta=ms2delta, by_col_e=if_merge_database_byenergy)
            db.normalize_peaks()

        db_frame = get_database_rankmz(db, maxrank=maxrank)
        db_frame.to_pickle(os.path.join(output_path, db.name + '_db_frame'+'.pklbz2'), compression='bz2', protocol=4)
