import sys
from sage.all import *
import cProfile
import random
import os
import math
import numpy
from scipy.stats import poisson


'''
This is the code that was run to load and process the data. There are many places where things can be tidied up, but this is intended to substantiate the claims in the paper, so I'm giving exactly what I ran.
'''


DEFAULT_DATA_PATH = ''
HDD_PATH_1 = '' # Fill in
HDD_PATH_2 = '' # Fill in


def main():
    pass


###########################
### Reading eigenvalues ###
###########################


def get_dirs_external_hdd(hdd_path=HDD_PATH_1):
    dir_names = ['10e5_p1009/', '300k_p1019/', '500k_p1013/',  '700k_p1019/',  '900k_p1019/', '10e5_p1013/', '300k_remaining/', '600k_p1009/', '700k_p1021/', '10e5-redo_p1009/', '300k_rerun/', '600k_p1013/', '800k_p1009/', '10e5-redo_p1013/', '400k_p1009/', '600k_p1019/', '800k_p1013/', 'redo400_p1009/', '300k_p1009/', '400k_p1013/', '700k_p1009/', '900k_p1009/', 'redo400_p1013/', '300k_p1013/', '500k_p1009/', '700k_p1013/', '900k_p1013/', 'redo400_p1019/']
    dir_names.sort()
    dir_paths = [hdd_path + name for name in dir_names]
    return dir_paths


def get_dirs_external_hdd_1M_2M(hdd_path=HDD_PATH_2):
    dir_names = ['11e5/', '12e5/', '13e5/', '14e5/', '15e5/', '16e5/', '17e5/', '18e5/', '19e5/', '20e5/']
    dir_names.sort()
    dir_paths = [hdd_path + name for name in dir_names]
    return dir_paths


def read_eigs_from_dirs(dir_list, num_eigs=300, levels=None, degs=None, run_remove_dupes=True, remove_cond=None, verbose=True, check_level_bounds_from_filename=False):
    hev = {}
    for directory in dir_list:
        filename_list = os.listdir(directory)
        if verbose:
            print(directory)
        for filename in sorted(filename_list):
            check_file = True
            if check_level_bounds_from_filename and (levels is not None):
                if 'heckeeigenvalues' in filename:
                    min_level_str = filename.split('_')[-2]
                    max_level_str = filename.split('_')[-1]
                    try:
                        min_level = int(min_level_str)
                        max_level = int(max_level_str)
                        any_levels_in_range = False
                        for level in levels:
                            if (min_level <= level) and (level <= max_level):
                                any_levels_in_range = True
                                break
                        if not any_levels_in_range:
                            check_file = False
                    except ValueError as e: # can't convert to int probably
                        print e
                        continue
            if check_file:
                new_hev = read_eigs(filename, num_eigs=num_eigs, levels=levels, degs=degs, data_path=directory)
                if run_remove_dupes:
                    new_hev = remove_dupes(new_hev)
                merge_hevs(hev, new_hev, verbose=verbose)
    if remove_cond:
        hev = get_hev_no_cond(hev, run_remove_dupes=run_remove_dupes)
    elif remove_cond is None:
        missing_cond_data = False
        for p in hev:
            k = hev[p].keys()
            if None in k:
                missing_cond_data = True
                break
        if missing_cond_data:
            hev = get_hev_no_cond(hev, run_remove_dupes=run_remove_dupes)
    return hev


def read_eigs(filename, num_eigs=None, levels=None, degs=None, data_path=DEFAULT_DATA_PATH):
    if 'heckeeigenvalues' not in filename:
        return {}
    f = open(data_path + filename,'r')
    lines = f.readlines()
    f.close()
    hev = {}
    z_minpoly = None
    conductor = None
    for line in lines:
        line = line.strip('\n')
        split_line = line.split(': ')
        
        if split_line[0] == 'p':
            p = int(split_line[1])
            if p not in hev:
                hev[p] = {}
            z_minpoly = None
            conductor = None
        
        if split_line[0] == 'conductor':
            conductor = int(split_line[1])
            if conductor not in hev[p]:
                hev[p][conductor] = []

        if split_line[0] == 'sign':
            conductor = int(split_line[1])
            if conductor not in hev[p]:
                hev[p][conductor] = []
        
        if split_line[0] == 'z.minpoly()':
            z_minpoly = ZZ['x'](sage_eval(split_line[1]))
        
        tmp_txt = 'Maximal Order in Number Field in z with defining polynomial '
        if (len(split_line[0]) > len(tmp_txt)) and (split_line[0][:len(tmp_txt)] == tmp_txt):
            z_minpoly = sage_eval(split_line[0][len(tmp_txt):], locals={'x':ZZ['x'].gen()})
        
        if (len(split_line[0]) > 0) and (split_line[0][0] == '['):
            level_ok = (levels is None) or (p in levels)
            deg_ok = (degs is None) or ((1 in degs) and (z_minpoly is None)) or ((z_minpoly is not None) and (z_minpoly.degree() in degs))
            if level_ok and deg_ok:
                if z_minpoly is not None:
                    K = NumberField(z_minpoly, 'z')
                    z = K.gen()
                    eigs = sage_eval(split_line[0], locals={'z':z})
                    if num_eigs is not None:
                        eigs = eigs[:num_eigs]
                    eigs = [K(tmp) for tmp in eigs]
                    z_minpoly = None
                else:
                    eigs = sage_eval(split_line[0])
                    if num_eigs is not None:
                        eigs = eigs[:num_eigs]
                    eigs = [QQ(tmp) for tmp in eigs]
                if conductor not in hev[p]:
                    hev[p][conductor] = []
                hev[p][conductor].append(eigs)

    if levels is not None:
        to_del = [p for p in hev if p not in levels]
        for p in to_del:
            del hev[p]
    
    return hev


def get_hev_no_cond(hev, run_remove_dupes=True):
    new_hev = {}
    for p in hev:
        new_hev[p] = []
        for conductor in hev[p]:
            new_hev[p] += hev[p][conductor]
    if run_remove_dupes:
        new_hev = remove_dupes_no_cond(new_hev)
    return new_hev


def remove_dupes(hev):
    # takes output from read_eigs
    new_hev = {}
    for p in hev:
        for conductor in hev[p]:
            if len(hev[p][conductor]) > 1:
                new_eigs_list = []
                traces_list = [[tmp.trace() for tmp in evl] for evl in hev[p][conductor]]
                for i,eigs in enumerate(hev[p][conductor]):
                    is_dupe = traces_list[i] in traces_list[i+1:]
                    if is_dupe:
                        print('Dupe at p = '+str(p))
                    else:
                        new_eigs_list.append(eigs)
                if p not in new_hev:
                    new_hev[p] = {}
                new_hev[p][conductor] = new_eigs_list
            else:
                if p not in new_hev:
                    new_hev[p] = {}
                new_hev[p][conductor] = hev[p][conductor]
    return new_hev


def remove_dupes_no_cond(hev):
    new_hev = {}
    for p in hev:
        new_eigs_list = []
        if len(hev[p]) > 1:
            traces_list = [[tmp.trace() for tmp in evl] for evl in hev[p]]
            for i,eigs in enumerate(hev[p]):
                is_dupe = traces_list[i] in traces_list[i+1:]
                if is_dupe:
                    print('Dupe at p = '+str(p))
                else:
                    new_eigs_list.append(eigs)
            new_hev[p] = new_eigs_list
        else:
            new_hev[p] = hev[p]
    return new_hev


def remove_shorter_dupes(f_list, verbose=False):
    if len(f_list) == 1:
        return f_list
    new_f_list = []
    traces_list = [[an.trace() for an in f] for f in f_list]
    for i, tl in enumerate(traces_list):
        good = True
        for j, other_tl in enumerate(traces_list):
            if i != j:
                len_tl = len(tl)
                len_other = len(other_tl)
                if len_tl < len_other:
                    is_sublist = True
                    for k in range(len_tl):
                        trc = tl[k]
                        other_trc = other_tl[k]
                        if trc != other_trc:
                            is_sublist = False
                            break
                    if is_sublist:
                        good = False
                        break
        if good:
            new_f_list.append(f_list[i])
        else:
            if verbose:
                print 'Removed shorter dupe'
    return new_f_list


def merge_hevs(old_hev, new_hev, verbose=True):
    for p in new_hev:
        if p not in old_hev:
            old_hev[p] = new_hev[p]
        else:
            for conductor in new_hev[p]:
                if conductor not in old_hev[p]:
                    old_hev[p][conductor] = new_hev[p][conductor]
                else:
                    old_traces = [[tmp.trace() for tmp in evl] for evl in old_hev[p][conductor]]
                    new_traces = [[tmp.trace() for tmp in evl] for evl in new_hev[p][conductor]]
                    
                    all_new_contained_in_old = True
                    for tl in new_traces:
                        if tl not in old_traces:
                            all_new_contained_in_old = False
                            break

                    if not all_new_contained_in_old:
                        all_old_contained_in_new = True
                        for tl in old_traces:
                            if tl not in new_traces:
                                all_old_contained_in_new = False
                                break
                        
                        if all_old_contained_in_new:
                            old_hev[p][conductor] = new_hev[p][conductor]
                        else:
                            if verbose:
                                print('not (all_old_contained_in_new or all_new_contained_in_old) at p = '+str(p))
                            for i,evl in enumerate(new_hev[p][conductor]):
                                tl = new_traces[i]
                                if tl not in old_traces:
                                    old_hev[p][conductor].append(evl)



#############################
### Reading deg sign disc ###
#############################


def read_level_to_deg_sign_disc(filename='level_to_deg_sign_disc'):
    #level_to_deg_sign_disc = load('level_to_deg_sign_disc.sobj')
    #return level_to_deg_sign_disc
    level_to_deg_sign_disc = {}
    f = open(filename,'r')
    for line in f:
        line = line.replace('\n','')
        split_line = line.split(':')
        level = int(split_line[0])
        deg_list = str_of_list_to_list(split_line[1], int)
        sign_list = str_of_list_to_list(split_line[2], int)
        disc_str = split_line[3]
        disc_str = disc_str.replace('?','-1')
        disc_list = str_of_list_to_list(disc_str, int)
        disc_list = disc_list[:-2] + [None,None]
        level_to_deg_sign_disc[level] = [(deg_list[i], sign_list[i], disc_list[i]) for i in range(len(deg_list))]
    f.close()
    return level_to_deg_sign_disc


def read_level_data(filename='level_data'):
    # Not very useful; use read_level_to_deg_sign_disc instead
    level_data = {}
    f = open(filename,'r')
    for line in f:
        line = line.replace('\n','')
        split_line = line.split(':')
        N = int(split_line[0])
        small_dims = str_of_list_to_list(split_line[1], int)
        discs = str_of_list_to_list(split_line[2], int)
        level_data[N] = [(small_dims[i],discs[i]) for i in range(len(small_dims))]
    f.close()
    return level_data


def read_sign_data(filename='atkin_lehner_data'):
    # Not very useful; use read_level_to_deg_sign_disc instead
    #al_data = load('al_data.sobj')
    f = open(filename,'r')
    for line in f:
        al_data = sage_eval(line) # takes about 10 seconds
        break
    f.close()
    return al_data


def get_setzer_neumann_levels():
    '''
    sage: SN_levels = []
    sage: for p in primes(10**4, 2*10**6):
    ....:     if sqrt(p - 64) in ZZ:
    ....:         SN_levels.append(p)
    ....:         
    sage: len(SN_levels)
    138
    sage: print(SN_levels)
    [15193, 16193, 18289, 21089, 21673, 23473, 26633, 27953, 30689, 31393, 37313, 38873, 42089, 45433, 51593, 59113, 60089, 65089, 69233, 70289, 71353, 75689, 80153, 91873, 93089, 99289, 104393, 106993, 112289, 119089, 120473, 124673, 134753, 140689, 142193, 156089, 162473, 164089, 165713, 170633, 191033, 199873, 214433, 216289, 225689, 227593, 233353, 235289, 237233, 247073, 267353, 277793, 284153, 286289, 288433, 305873, 310313, 319289, 321553, 330689, 332993, 368513, 375833, 378289, 380753, 400753, 403289, 416089, 444953, 466553, 508433, 511289, 543233, 552113, 585289, 597593, 600689, 603793, 613153, 616289, 628913, 632089, 664289, 667553, 700633, 710713, 727673, 744833, 751753, 786833, 822713, 833633, 859393, 889313, 912089, 931289, 950689, 966353, 986113, 994073, 1014113, 1034353, 1075433, 1092089, 1113089, 1155689, 1159993, 1172953, 1181633, 1199089, 1221089, 1238833, 1270193, 1283753, 1315673, 1361953, 1385393, 1399553, 1404289, 1409033, 1447273, 1471433, 1481153, 1530233, 1570073, 1646153, 1677089, 1697873, 1703089, 1724033, 1776953, 1782289, 1787633, 1814473, 1841513, 1868753, 1923833, 1979713]
    '''
    to_ret = [15193, 16193, 18289, 21089, 21673, 23473, 26633, 27953, 30689, 31393, 37313, 38873, 42089, 45433, 51593, 59113, 60089, 65089, 69233, 70289, 71353, 75689, 80153, 91873, 93089, 99289, 104393, 106993, 112289, 119089, 120473, 124673, 134753, 140689, 142193, 156089, 162473, 164089, 165713, 170633, 191033, 199873, 214433, 216289, 225689, 227593, 233353, 235289, 237233, 247073, 267353, 277793, 284153, 286289, 288433, 305873, 310313, 319289, 321553, 330689, 332993, 368513, 375833, 378289, 380753, 400753, 403289, 416089, 444953, 466553, 508433, 511289, 543233, 552113, 585289, 597593, 600689, 603793, 613153, 616289, 628913, 632089, 664289, 667553, 700633, 710713, 727673, 744833, 751753, 786833, 822713, 833633, 859393, 889313, 912089, 931289, 950689, 966353, 986113, 994073, 1014113, 1034353, 1075433, 1092089, 1113089, 1155689, 1159993, 1172953, 1181633, 1199089, 1221089, 1238833, 1270193, 1283753, 1315673, 1361953, 1385393, 1399553, 1404289, 1409033, 1447273, 1471433, 1481153, 1530233, 1570073, 1646153, 1677089, 1697873, 1703089, 1724033, 1776953, 1782289, 1787633, 1814473, 1841513, 1868753, 1923833, 1979713]
    return to_ret


def get_level_to_deg_sign_disc_no_setzer_neumann():
    ldsd = read_level_to_deg_sign_disc()
    new_ldsd = {}
    for level in ldsd:
        sn_ok = not (sqrt(level - 64) in ZZ)
        new_v = []
        for (tmp_deg, tmp_sign, tmp_disc) in ldsd[level]:
            if sn_ok or (not (tmp_deg, tmp_sign) == (1,-1)):
                new_v.append((tmp_deg, tmp_sign, tmp_disc))
            else:
                sn_ok = True
        new_ldsd[level] = new_v
    return new_ldsd


def hev_to_field_to_levels(hev):
    # hev output of read_eigs_from_dirs, e.g.
    field_to_levels = {}
    for p in sorted(hev.keys()):
        if isinstance(hev[p], dict):
            for conductor in hev[p]:
                for eigs in hev[p][conductor]:
                    K = get_K(eigs)
                    if K not in field_to_levels:
                        field_to_levels[K] = []
                    field_to_levels[K].append(p)
        else:
            for eigs in hev[p]:
                K = get_K(eigs)
                if K not in field_to_levels:
                    field_to_levels[K] = []
                field_to_levels[K].append(p)
    return field_to_levels


def get_space_sizes(al_signs):
    '''
    al_signs dict {level:[(dim1, sign1), ...]}
    '''
    space_sizes = {}
    for N in al_signs:
        plus_size = 0
        minus_size = 0
        for (deg, sgn) in al_signs[N]:
            if sgn == 1:
                plus_size += deg
            else:
                minus_size += deg
        space_sizes[N] = (plus_size, minus_size)
    return space_sizes



############################
### Making deg sign disc ###
############################


def make_level_data_with_signs(level_data, al_data, verbose=True):
    level_data_with_signs = {} # {level:[(deg, sign, disc)]}
    bad_levels = set()
    levels_done_by_hand = [55949, 62297, 85829, 131627, 367853] # Done in heckeeigenvalues_2022-12-14_55949_367853
    for level in primes(10**4, 2*10**6):
        if level not in levels_done_by_hand:
            new_ld = []
            if len(al_data[level]) > 2:
                deg_to_discs = {}
                for (deg,disc) in level_data[level]:
                    if deg not in deg_to_discs:
                        deg_to_discs[deg] = []
                    deg_to_discs[deg].append(disc)
                deg_to_signs = {}
                for (deg,sgn) in al_data[level][:-2]:
                    if deg not in deg_to_signs:
                        deg_to_signs[deg] = []
                    deg_to_signs[deg].append(sgn)

                for deg in deg_to_discs:
                    num_discs = len(set(deg_to_discs[deg]))
                    num_signs = len(set(deg_to_signs[deg]))
                    if num_discs == 1:
                        disc = deg_to_discs[deg][0]
                        for (tmp_deg,tmp_sgn) in al_data[level]:
                            if tmp_deg == deg:
                                new_datum = (deg, tmp_sgn, disc)
                                new_ld.append(new_datum)
                    elif num_signs == 1:
                        sgn = deg_to_signs[deg][0]
                        for (tmp_deg,tmp_disc) in level_data[level]:
                            if tmp_deg == deg:
                                new_datum = (deg, sgn, tmp_disc)
                                new_ld.append(new_datum)
                    else:
                        bad_levels.add(level)
                        if verbose:
                            message = '### Can\'t match signs and discs at level '+str(level)+', degree '+str(deg)+'.'
                            print message
                            print 'al_data:', al_data[level]
                            print 'level_data:', level_data[level]

            for (deg, sgn) in al_data[level][-2:]:
                new_datum = (deg, sgn, None)
                new_ld.append(new_datum)
        else:
            if level == 55949:
                new_ld = [(1, 1, 1), (2, -1, 5), (2, 1, 5), (2, 1, 13), (2247, 1, None), (2408, -1, None)]
            elif level == 62297:
                new_ld = [(2, 1, 5), (2, -1, 8), (2546, 1, None), (2641, -1, None)]
            elif level == 85829:
                new_ld = [(2, 1, 5), (2, -1, 8), (3480, 1, None), (3668, -1, None)]
            elif level == 131627:
                new_ld = [(1, 1, 1), (1, -1, 1), (1, 1, 1), (2, 1, 5), (2, -1, 8), (5370, 1, None), (5592, -1, None)]
            elif level == 367853:
                new_ld = [(2, 1, 5), (2, -1, 8), (15177, 1, None), (15473, -1, None)]
            else:
                raise NotImplementedError
        
        new_ld.sort()
        level_data_with_signs[level] = new_ld
    to_ret = (level_data_with_signs, bad_levels)
    return to_ret


def level_with_sign_disc_to_file(level_to_sign_disc, filename):
    to_write = ''
    for level in sorted(level_to_sign_disc.keys()):
        deg_list = []
        sign_list = []
        disc_list = []
        for (deg,sgn,disc) in level_to_sign_disc[level]:
            deg_list.append(deg)
            sign_list.append(sgn)
            disc_list.append(disc)
        level_str = str(level)
        deg_str = str(deg_list)
        sign_str = str(sign_list)
        disc_str = str(disc_list)
        disc_str = disc_str.replace('None','?')
        to_write += level_str + ':' + deg_str + ':' + sign_str + ':' + disc_str + '\n'
    to_write = to_write[:-1]
    to_write = to_write.replace(' ','')
    f = open(filename,'w')
    f.write(to_write)
    f.close()
    return


def charpoly_file_to_dims(filename, data_path=DEFAULT_DATA_PATH):
    if 'charpolys' not in filename:
        return {}
    f = open(data_path + filename,'r')
    lines = f.readlines()
    f.close()
    dims = {}
    for line in lines:
        line = line.strip('\n')
        split_line = line.split(': ')
        
        if split_line[0] == 'p':
            p = int(split_line[1])
            if p not in dims:
                dims[p] = set()

        if split_line[0] == 'minpoly':
            minpoly_deg = split_line[1].count(',') # number of commas is length of list -1, and that's degree of poly

        if split_line[0] == 'repeatedfactors':
            repeatedfactors_deg = -1 + split_line[1].count('[')
            total_deg = minpoly_deg + repeatedfactors_deg
            dims[p].add(total_deg)

    for p,v in dims.items():
        new_v = list(v)
        new_v.sort()
        try:
            new_v[1] -= 1 # remove x-3 Eisenstein factor
        except IndexError:
            print 'strange dims at level ', p, ':', new_v
        new_v = tuple(new_v)
        dims[p] = new_v
    
    return dims


def dirs_to_dims(dir_list, verbose=True):
    dims = {}
    for directory in dir_list:
        filename_list = os.listdir(directory)
        if verbose:
            print(directory)
        for filename in sorted(filename_list):
            new_dims = charpoly_file_to_dims(filename, data_path=directory)
            for N in new_dims:
                if (N not in dims) or (len(dims[N]) < len(new_dims[N])):
                    dims[N] = new_dims[N]
    return dims


def magma_file_to_signs(filename):
    signs_dict = {}
    f = open(filename,'r')
    for line in f:
        # N:k:i:t:D:T:A:F:C:E:cm:tw:pra:zr:mm:h:X:sd:eap
        split_line = line.split(':')
        N = int(split_line[0])
        D = str_of_list_to_list(split_line[4], int)
        A_str = split_line[6]
        A_str_split = A_str.split('>')
        signs = []
        for tmpstr in A_str_split[:-1]:
            signstr = tmpstr[-2:]
            signstr = signstr.replace(',','')
            signs.append(int(signstr))
        deg_to_sign = [(D[i],signs[i]) for i in range(len(D))]
        signs_dict[N] = deg_to_sign
    return signs_dict


def magma_dir_to_signs(data_path):
    dirlist = os.listdir(data_path)
    signs_dict = {}
    for filename in dirlist:
        print filename
        new_signs_dict = magma_file_to_signs(data_path + filename)
        signs_dict.update(new_signs_dict)
    return signs_dict



####################
### Data to file ###
####################


def make_datapoints_files():
    tmp = datapoints_uptoX_by_deg(degs=[1,2,3,4,5,6], to_file=True)
    tmp = datapoints_uptoX_deg2_by_disc(disclist_deg2=[5,8,12,13,17,21], to_file=True)
    tmp = datapoints_uptoX_deg3_by_disc(disclist_deg3=[49,81,148,169,229,257,321], to_file=True)
    tmp = datapoints_delta(to_file=True)
    tmp = datapoints_deg3_C2_histo(to_file=True)
    tmp = datapoints_deg3_C_23(to_file=True)
    tmp = datapoints_deg2_C23(to_file=True) # This takes 6 GB of RAM!
    tmp = datapoints_poisson_logprobs(to_file=True)
    return


def datapoints_uptoX_by_deg(degs=[1,2,3,4,5,6], to_file=False):
    counts = make_counts_by_deg(degs=degs)
    if to_file:
        for d in degs:
            filename = 'uptoX_deg_'+str(d)
            write_str = ''
            for (level, c) in counts[d]:
                write_str += str(level)+','+str(c)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()
    return counts


def datapoints_uptoX_deg2_by_disc(disclist_deg2=[5,8,12,13,17,21], to_file=False):
    level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    starting_counts_deg2 = {5:158, 8:37, 12:1, 13:13, 17:0, 21:1}
    counts_deg2 = {disc: make_counts_by_deg(degs=[2], disc_list=[disc], starting_counts={2:starting_counts_deg2[disc]})[2] for disc in disclist_deg2}
    if to_file:
        for disc in disclist_deg2:
            filename = 'uptoX_deg_2_disc_'+str(disc)
            write_str = ''
            for (level, c) in counts_deg2[disc]:
                write_str += str(level)+','+str(c)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()
    return counts_deg2


def datapoints_uptoX_deg3_by_disc(disclist_deg3=[49,81,148,169,229,257,321], to_file=False):
    level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    starting_counts_deg3 = {49:34, 81:3, 148:12, 169:2, 229:8, 257:9, 321:2}
    counts_deg3 = {disc: make_counts_by_deg(degs=[3], disc_list=[disc], starting_counts={3:starting_counts_deg3[disc]})[3] for disc in disclist_deg3}
    if to_file:
        for disc in disclist_deg3:
            filename = 'uptoX_deg_3_disc_'+str(disc)
            write_str = ''
            for (level, c) in counts_deg3[disc]:
                write_str += str(level)+','+str(c)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()
    return counts_deg3


def datapoints_delta(to_file=False):
    level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    delta_1 = get_delta_by_degdisclist(level_to_deg_sign_disc, 1, disc_list=None)
    delta_1_nosn = get_delta_by_degdisclist(level_to_deg_sign_disc, 1, disc_list=None, omit_SN=True)
    delta_2 = get_delta_by_degdisclist(level_to_deg_sign_disc, 2, disc_list=None)
    delta_3 = get_delta_by_degdisclist(level_to_deg_sign_disc, 3, disc_list=None)
    
    delta_1_list = sorted(delta_1.items())
    delta_1_nosn_list = sorted(delta_1_nosn.items())
    delta_2_list = sorted(delta_2.items())
    delta_3_list = sorted(delta_3.items())
    
    delta_1_list = [(k,int(round(v))) for k,v in delta_1_list]
    delta_1_nosn_list = [(k,int(round(v))) for k,v in delta_1_nosn_list]
    delta_2_list = [(k,int(round(v))) for k,v in delta_2_list]
    delta_3_list = [(k,int(round(v))) for k,v in delta_3_list]
    
    to_ret = [delta_1_list, delta_1_nosn_list, delta_2_list, delta_3_list]
    
    if to_file:
        for i,deltalist in enumerate(to_ret):
            filename = 'heads_minus_tails_deg_'+str(['1','1_nosn','2','3'][i])
            write_str = ''
            for (level, c) in deltalist:
                write_str += str(level)+','+str(c)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()

    return to_ret


def datapoints_deg3_C2_histo(to_file=False):
    to_ret = {
        49: [2, 2, 2, 1, 0, 1, 1, 1, 0, 1, 0, 3, 0, 5, 0, 1, 1, 2, 2, 3, 3, 0, 1, 2, 1, 0, 1, 1, 1, 2, 2, 2, 0, 0, 1, 3, 1, 3, 2, 2, 2, 1, 3, 1, 2, 2, 0, 3, 2, 0, 1, 3, 2, 2, 2, 4, 0, 0, 4, 1, 4, 1, 1, 0, 4, 2, 0, 1, 3, 2, 2, 1, 0, 2, 3, 2, 1, 1, 2, 2, 1, 5, 1, 3, 3, 4, 1, 1, 2, 2, 0, 2, 1, 1, 2, 0, 1, 3, 0, 2, 0, 4, 0, 1, 3, 1, 3, 7, 3, 1, 1, 1, 0, 0, 2, 1, 1, 2, 0, 1],
        81: [4, 0, 4, 0, 1, 1, 0, 1, 3, 2, 3, 3, 4],
        148: [3, 2, 5, 7, 3, 4],
        169: [1, 2, 3, 4, 0, 0, 2, 1, 7],
        229: [0, 2, 2, 4, 2, 3, 2, 4, 1, 1, 2, 3, 3, 4, 4, 3, 1, 2, 2, 1, 4],
        257: [16, 6, 17, 6, 4, 1, 15],
        321: [7],
    }
    if to_file:
        for disc,vallist in sorted(to_ret.items()):
            filename = 'coincidences_k_2_deg_3_disc_'+str(disc)+'_values'
            write_str = ''
            for val in vallist:
                write_str += str(val)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()

    return to_ret


def datapoints_deg3_C_23(to_file=False):
    '''
    pia_counts_3 = load('pia_counts_3.sobj')
    formlens_3 = load('formlens_3.sobj')
    counts_by_D = {2:{}, 3:{}}
    for level in pia_counts_3:
        for i in range(len(pia_counts_3[level])):
            D = formlens_3[level][i][1].number_field().disc()
            if D not in counts_by_D[2]:
                counts_by_D[2][D] = []
            if D not in counts_by_D[3]:
                counts_by_D[3][D] = []
            c2 = pia_counts_3[level][i].get(2,0)
            c3 = pia_counts_3[level][i].get(3,0)
            counts_by_D[2][D].append((level,c2))
            counts_by_D[3][D].append((level,c3))
    return counts_by_D
    '''
    pia_counts_3 = load('pia_counts_3.sobj')
    formlens_3 = load('formlens_3.sobj')
    form_dicts = (formlens_3, pia_counts_3)
    counts = form_dicts_to_pia_points(form_dicts, countvals=[2,3], by_disc=True, normalize=False)
    if to_file:
        for disc in sorted(counts.keys()):
            for c in sorted(counts[disc].keys()):
                filename = 'coincidences_k_'+str(c)+'_deg_3_disc_'+str(disc)
                write_str = ''
                for (level, c) in counts[disc][c]:
                    write_str += str(level)+','+str(c)+'\n'
                write_str = write_str[:-1]
                f = open(filename,'w')
                f.write(write_str)
                f.close()
    return counts
    

def get_c2_badtups():
    # Some of the early q-expansions weren't computed far enough. They've been recomputed, but I haven't removed the old ones yet.
    '''
    sage: for d in c22_points:
    ....:     print ''
    ....:     print d
    ....:     for k,v in c22_points[d]:
    ....:         maybedupes = [tup for tup in c22_points[d] if tup[0] == k]
    ....:         expected = len([tup for tup in ldsd[k] if tup[0]==2 and tup[2]==d])
    ....:         if len(maybedupes) != expected:
    ....:             print maybedupes, ldsd[k], expected==1
    ....:             badtups.append(maybedupes[0])
    ....:             
    
    5
    [(89591, 37), (89591, 42)] [(2, 1, 5), (3539, 1, None), (3925, -1, None)] True
    [(89591, 37), (89591, 42)] [(2, 1, 5), (3539, 1, None), (3925, -1, None)] True
    [(169943, 29), (169943, 90)] [(2, -1, 5), (6958, 1, None), (7202, -1, None)] True
    [(169943, 29), (169943, 90)] [(2, -1, 5), (6958, 1, None), (7202, -1, None)] True
    [(179327, 21), (179327, 81)] [(2, 1, 5), (7339, 1, None), (7603, -1, None)] True
    [(179327, 21), (179327, 81)] [(2, 1, 5), (7339, 1, None), (7603, -1, None)] True
    [(191803, 32), (191803, 97)] [(1, 1, 1), (2, -1, 5), (7912, 1, None), (8068, -1, None)] True
    [(191803, 32), (191803, 97)] [(1, 1, 1), (2, -1, 5), (7912, 1, None), (8068, -1, None)] True
    [(309931, 112), (309931, 115)] [(1, 1, 1), (2, -1, 5), (12772, 1, None), (13052, -1, None)] True
    [(309931, 112), (309931, 115)] [(1, 1, 1), (2, -1, 5), (12772, 1, None), (13052, -1, None)] True
    [(387187, 110), (387187, 125)] [(2, 1, 5), (16066, 1, None), (16197, -1, None)] True
    [(387187, 110), (387187, 125)] [(2, 1, 5), (16066, 1, None), (16197, -1, None)] True

    8
    
    12
    
    13
    
    17
    
    21
    '''
    badtups = [(89591, 37), (169943, 29), (179327, 21), (191803, 32), (309931, 112), (387187, 110)]
    return badtups


def datapoints_deg2_C23(to_file=False):
    # This takes 6 GB of RAM!
    dtu2 = load('deg2_formlens_pia_piaall_ratans_ratansall_10k_1M.sobj')
    dtu2_2 = load('deg2_formlens_pia_piaall_ratans_ratansall_1M_2M.sobj')
    for i in range(5):
        dtu2[i].update(dtu2_2[i])
    counts = form_dicts_to_pia_points(dtu2, countvals=[2,3]) # this will have the badtups
    badtups = get_c2_badtups()
    if to_file:
        for disc in sorted(counts.keys()):
            for c in sorted(counts[disc].keys()):
                filename = 'coincidences_k_'+str(c)+'_deg_2_disc_'+str(disc)
                write_str = ''
                for (level, cc) in counts[disc][c]:
                    if (level, cc) not in badtups:
                        write_str += str(level)+','+str(cc)+'\n'
                write_str = write_str[:-1]
                f = open(filename,'w')
                f.write(write_str)
                f.close()
    return counts


def datapoints_poisson_logprobs(to_file=False):
    lp1 = load('poisson_logprobs_deg1.sobj')
    lp5 = load('poisson_logprobs_deg2_disc5.sobj')
    lp8 = load('poisson_logprobs_deg2_disc8.sobj')
    lp49 = load('poisson_logprobs_deg3_disc49.sobj')
    to_ret = [lp1,lp5,lp8,lp49]
    if to_file:
        for i,lp in enumerate(to_ret):
            filename = 'loglikelihoods_deg_'+str([1,2,2,3][i])+'_disc_'+str([1,5,8,49][i])
            write_str = ''
            for ((a,b),p) in sorted(lp.items()):
                write_str += str(a)+','+str(b)+','+str(p)+'\n'
            write_str = write_str[:-1]
            f = open(filename,'w')
            f.write(write_str)
            f.close()
    return to_ret



############
### Misc ###
############


def str_of_list_to_list(str_list, entry_type):
    to_ret = []
    str_list = str_list.strip('[]')
    if str_list == '':
        return to_ret
    if ',' not in str_list:
        to_ret.append(entry_type(str_list))
        return to_ret
    str_entries = str_list.split(',')
    for entry in str_entries:
        to_ret.append(entry_type(entry))
    return to_ret


def get_K(eigs):
    K = QQ
    for a in eigs:
        if (a not in K): #or (a not in ZZ):
            K = a.parent()
            break
    if hasattr(K, 'number_field'):
        K = K.number_field()
    return K


def get_expected_fn_power(power=1):
    def expected_fn(plus_dim, minus_dim):
        return float(minus_dim)**power/(plus_dim**power + minus_dim**power)
    return expected_fn


def get_points_kwargs(index, label, size, small=False, marker=None, s=1, v=1, color_shift=0, zorder=0, extra=None):
    color = get_colour_list(s=s, v=v, shift=color_shift)[index]
    if small:
        zorder += 1
    faceted = small
    if faceted:
        markeredgecolor = 'black'
    else:
        markeredgecolor = None
    pt_kwargs = {'color':color, 'size':size, 'faceted':faceted, 'markeredgecolor':markeredgecolor, 'legend_label':label, 'legend_color':color, 'marker':marker, 'zorder':zorder}
    if extra is not None:
        pt_kwargs.update(extra)
    return pt_kwargs


def get_colour_list(s=1, v=1, shift=0):
    colour_list = [hue(0.58+shift, s=s, v=v), hue(0.70+shift, s=s, v=v), hue(0.8+shift, s=s, v=v), hue(0.0+shift, s=s, v=v), hue(0.12+shift, s=s, v=v), hue(0.22+shift, s=s, v=v), hue(0.46+shift, s=s, v=v)]
    #colour_list = [hue(0.58+shift, s=s, v=v), hue(0.70+shift, s=s, v=v), hue(0.8+shift, s=s, v=v), hue(0.0+shift, s=s, v=v), hue(0.12+shift, s=s, v=v), hue(0.26+shift, s=s, v=v), hue(0.46+shift, s=s, v=v)]
    #colour_list = [hue(0.50+shift, s=s, v=v), hue(0.62+shift, s=s, v=v), hue(0.76+shift, s=s, v=v), hue(0.0+shift, s=s, v=v), hue(0.12+shift, s=s, v=v), hue(0.26+shift, s=s, v=v), hue(0.46+shift, s=s, v=v)]
    return colour_list


# "The ratio of the Plancherel measures of these sets of integers is about 45:12:19:24."
'''
sage: dat = [(1, [x - 2, x - 1, x, x + 1, x + 2]),
....:  (13, [x^2 - x - 3, x^2 + x - 3]),
....:  (13^2,
....:   [x^3 - 2*x^2 - 3*x + 5,
....:    x^3 - x^2 - 4*x - 1,
....:    x^3 + x^2 - 4*x + 1,
....:    x^3 + 2*x^2 - 3*x - 5]),
....:  (13^5,
....:   [x^6 - 5*x^5 + 5*x^4 + 6*x^3 - 7*x^2 - 2*x + 1,
....:    x^6 - x^5 - 5*x^4 + 4*x^3 + 6*x^2 - 3*x - 1,
....:    x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1,
....:    x^6 + 5*x^5 + 5*x^4 - 6*x^3 - 7*x^2 + 2*x + 1])]
sage: dat = dict(dat)

sage: def get_thetas(a2poly):
....:     a2poly = RR['x'](a2poly)
....:     roots = []
....:     for r,m in a2poly.roots():
....:         for _ in range(m):
....:             roots.append(r)
....:     lmbdas = [r/sqrt(2.0) for r in roots]
....:     thetas = [arccos(l/2) for l in lmbdas]
....:     return thetas

sage: def plancherel(theta, p=2):
....:     f1 = 2/RR(pi)
....:     f2 = (p+1)*sin(theta)**2
....:     t1 = (p**(0.5) + p**(-0.5))**2
....:     t2 = 4*cos(theta)**2
....:     f3 = t1 - t2
....:     to_ret = f1*f2/f3
....:     return to_ret

sage: plancherel_vals = {k:[] for k in dat}
sage: for k,v in sorted(dat.items()):
....:     print ''
....:     print k
....:     for poly in v:
....:         #print ''
....:         print poly
....:         thetas = get_thetas(poly)
....:         plancherel_prod = 1
....:         for t in thetas:
....:             plancherel_prod *= plancherel(t)
....:         plancherel_prod = plancherel_prod**(6/poly.degree())
....:         print plancherel_prod
....:         plancherel_vals[k].append(plancherel_prod)
....:         

1
x - 2
0.00310590551667784
x - 1
0.00531736700405063
x
0.00584430918329191
x + 1
0.00531736700405063
x + 2
0.00310590551667784

13
x^2 - x - 3
0.00295735724438803
x^2 + x - 3
0.00295735724438803

169
x^3 - 2*x^2 - 3*x + 5
0.00307779085150749
x^3 - x^2 - 4*x - 1
0.00167503547029818
x^3 + x^2 - 4*x + 1
0.00167503547029818
x^3 + 2*x^2 - 3*x - 5
0.00307779085150749

371293
x^6 - 5*x^5 + 5*x^4 + 6*x^3 - 7*x^2 - 2*x + 1
0.00133479719338609
x^6 - x^5 - 5*x^4 + 4*x^3 + 6*x^2 - 3*x - 1
0.00461858124358234
x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1
0.00461858124358234
x^6 + 5*x^5 + 5*x^4 - 6*x^3 - 7*x^2 + 2*x + 1
0.00133479719338609

sage: sumvals = [sum(plancherel_vals[k]) for k in plancherel_vals]
sage: sumvals_n = [tmp/sum(sumvals) for tmp in sumvals]
sage: sumvals_n
[0.453653966578209, 0.118251770622380, 0.190044719514594, 0.238049543284817]
sage: 
'''

##############
### Counts ###
##############


def get_counts_up_to_X_by_degdisc(level_to_deg_sign_disc, deg, disc, expected_fn=None, omit_SN=False, min_level=None, size=2):
    if expected_fn is None:
        def expected_fn(tmp1,tmp2):
            return 0.5
    
    counts = [(0,0)]
    for level in sorted(level_to_deg_sign_disc.keys()):
        if (min_level is None) or (level > min_level):
            if omit_SN and (deg == 1) and (sqrt(level - 64) in ZZ):
                SN_ok = False
            else:
                SN_ok = True
            for tmp_deg, tmp_sign, tmp_disc in level_to_deg_sign_disc[level]:
                if SN_ok or ((tmp_deg, tmp_sign) != (1,-1)):
                    if (tmp_deg == deg) and (tmp_disc == disc):
                        plus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == -1])
                        counts.append((level, counts[-1][1] + expected_fn(plus_dim, minus_dim) * 2))
                else:
                    SN_ok = True
    to_ret = counts
    return to_ret


def get_plusminus_counts_up_to_X(al_data, deg, with_plot=True, expected_fn=None, omit_SN=False, min_level=None, size=2):
    if expected_fn is None:
        def expected_fn(tmp1,tmp2):
            return 0.5
    
    plus_counts = [(0,0)]
    minus_counts = [(0,0)]
    for level in sorted(al_data.keys()):
        if (min_level is None) or (level > min_level):
            if omit_SN and (deg == 1) and (sqrt(level - 64) in ZZ):
                SN_ok = False
            else:
                SN_ok = True
            for d, s in al_data[level]:
                if SN_ok or ((d,s) != (1,-1)):
                    if d == deg:
                        plus_dim = sum([tmp[0] for tmp in al_data[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in al_data[level] if tmp[1] == -1])
                        if s == 1:
                            plus_counts.append((level, plus_counts[-1][1] + expected_fn(plus_dim, minus_dim) * 2))
                        else:
                            minus_counts.append((level, minus_counts[-1][1] + (1 - expected_fn(plus_dim, minus_dim)) * 2))
                else:
                    SN_ok = True
    to_ret = [plus_counts, minus_counts]
    if with_plot:
        plusplt = scatter_plot(plus_counts, markersize=size, facecolor='red', edgecolor='red')
        minusplt = scatter_plot(minus_counts, markersize=size, facecolor='blue', edgecolor='blue')
        plt = plusplt + minusplt
        show(plt)
        to_ret += [plt]
    to_ret = tuple(to_ret)
    return to_ret


def get_plusminus_counts_up_to_X_by_degdisc(level_to_deg_sign_disc, deg, disc, with_plot=True, expected_fn=None, omit_SN=False, min_level=None, size=2):
    if expected_fn is None:
        def expected_fn(tmp1,tmp2):
            return 0.5
    
    plus_counts = [(0,0)]
    minus_counts = [(0,0)]
    for level in sorted(level_to_deg_sign_disc.keys()):
        if (min_level is None) or (level > min_level):
            if omit_SN and (deg == 1) and (sqrt(level - 64) in ZZ):
                SN_ok = False
            else:
                SN_ok = True
            for tmp_deg, tmp_sign, tmp_disc in level_to_deg_sign_disc[level]:
                if SN_ok or ((tmp_deg, tmp_sign) != (1,-1)):
                    if (tmp_deg == deg) and (tmp_disc == disc):
                        plus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == -1])
                        if tmp_sign == 1:
                            plus_counts.append((level, plus_counts[-1][1] + expected_fn(plus_dim, minus_dim) * 2))
                        else:
                            minus_counts.append((level, minus_counts[-1][1] + (1 - expected_fn(plus_dim, minus_dim)) * 2))
                else:
                    SN_ok = True
    to_ret = [plus_counts, minus_counts]
    if with_plot:
        plusplt = points(plus_counts, size=size, color='red')
        minusplt = points(minus_counts, size=size, color='blue')
        plt = plusplt + minusplt
        show(plt)
        to_ret += [plt]
    to_ret = tuple(to_ret)
    return to_ret


def get_delta(al_data, degs, omit_SN=False, min_level=10**4, expected_fn=None):
    if expected_fn is None:
        def expected_fn(tmp1,tmp2):
            return 0.5
    
    count_diffs = {d:{} for d in degs}
    diff_tally = {d:0 for d in degs}
    for level,v in sorted(al_data.items()):
        if level > min_level:
            sn_ok = (not omit_SN) or (not (sqrt(level - 64) in ZZ))
            for deg,sgn in v:
                if deg in count_diffs:
                    if sn_ok or (not (deg,sgn) == (1,-1)):
                        plus_dim = sum([tmp[0] for tmp in al_data[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in al_data[level] if tmp[1] == -1])
                        if sgn == 1:
                            diff_tally[deg] += expected_fn(plus_dim, minus_dim) * 2
                        else:
                            diff_tally[deg] -= (1 - expected_fn(plus_dim, minus_dim)) * 2
                        count_diffs[deg][level] = diff_tally[deg]
                    else:
                        sn_ok = True
    return count_diffs


def make_deg2_up_to_X_graphs_by_disc(level_to_deg_sign_disc=None):
    if level_to_deg_sign_disc is None:
        level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    disclist_deg2 = [5,8,12,13,17,21]
    counts_deg2 = {disc: get_counts_up_to_X_by_degdisc(level_to_deg_sign_disc, 2, disc) for disc in disclist_deg2}
    pltlist = []
    for i,disc in enumerate(disclist_deg2):
        if i==0:
            size = 4
        elif i==1:
            size = 5
        else:
            size = 30
        label = str(disc)
        small = len(counts_deg2[disc]) < 20
        if i == 0:
            color_shift = -0.05
        else:
            color_shift = 0
        plt_kwargs = get_points_kwargs(i, label, size, small, extra={'ymin':0, 'xmax':2e6, 'dpi':800}, color_shift=color_shift)
        plt = points(counts_deg2[disc][1:], **plt_kwargs)
        if i < 2:
            plt.set_legend_options(font_size='large', markerscale=4.3, loc=(0.03,0.66))
        else:
            plt.set_legend_options(font_size='large', markerscale=1.0, loc=(0.025,0.80))
        pltlist.append(plt)
    plt1 = sum(pltlist)
    plt1.set_legend_options(loc=(0.03,0.73))
    to_ret = (plt1, sum(pltlist[2:]))
    return to_ret


def make_deg3_up_to_X_graphs_by_disc(level_to_deg_sign_disc=None, noise=0, dpi=1000):
    if level_to_deg_sign_disc is None:
        level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    disclist_deg3 = [49,81,148,169,229,257,321]
    counts_deg3 = {disc: get_counts_up_to_X_by_degdisc(level_to_deg_sign_disc, 3, disc) for disc in disclist_deg3}
    if noise != 0:
        counts_deg3 = {disc: [(tmp[0], tmp[1] + noise*random.random()) for tmp in counts_deg3[disc]] for disc in counts_deg3}
    pltlist = []
    pltlist2 = []
    for i,disc in enumerate(disclist_deg3):
        size = 20
        label = str(disc)
        small = len(counts_deg3[disc]) < 10
        zorder = 10**10 - len(counts_deg3[disc])
        extra = {'ymin':0, 'xmax':2e6, 'dpi':dpi}
        extra2 = dict(extra)
        extra2.update({'scale':'semilogx'})
        plt_kwargs = get_points_kwargs(i, label, size, small, zorder=zorder, extra=extra)
        plt_kwargs2 = get_points_kwargs(i, label, size, small, zorder=zorder, extra=extra2)
        plt = points(counts_deg3[disc][1:], **plt_kwargs)
        plt2 = points(counts_deg3[disc][1:], **plt_kwargs2)
        plt.set_legend_options(font_size='large', markerscale=1.0, loc=(0.03, 0.69))
        plt2.set_legend_options(font_size='large', markerscale=1.0, loc=(0.00, 0.72))
        pltlist.append(plt)
        if i > 0:
            pltlist2.append(plt2)
    plt1 = sum(pltlist)
    plt2 = sum(pltlist2)
    to_ret = (plt1, plt2)
    return to_ret



#########################
### Fitting to counts ###
#########################


def get_fits(degs, spacing=1, show_plots=True, loglog=False):
    counts_by_deg = make_counts_by_deg(degs=degs, loglog=loglog)
    model = get_model(loglog=loglog)
    fits = {}
    plots = {}
    for d in degs:
        fits[d] = get_fit_from_counts(counts_by_deg, d, model, spacing=spacing)
        if loglog:
            scale=None
        else:
            scale='loglog'
        plots[d] = make_plot(counts_by_deg, d, model, fits[d], show_plot=False, scale=scale)
    if show_plots:
        for d in degs:
            print d,':',fits[d]
        show(sum(flatten(plots.values())))
    to_ret = (plots, fits)
    return to_ret


def get_fits_by_disc(deg, disc_list, spacing=100, show_plots=True, loglog=True):
    level_to_deg_sign_disc = read_level_to_deg_sign_disc()
    counts = {disc: get_counts_up_to_X_by_degdisc(level_to_deg_sign_disc, deg, disc) for disc in disc_list}
    starting_counts = {1:329, 5:158, 8:37, 12:1, 13:13, 17:0, 21:1, 49:34, 229:8, 148:12, 81:3, 257:9, 169:2, 321:2, 725:16, 1957:4, 2777:3, 8768:0, 70601:2, 11**4:0, 13**5:0}
    new_counts = {}
    for disc,tuplist in counts.items():
        new_tuplist = [(tmp[0],tmp[1] + starting_counts[disc]) for tmp in tuplist[1:]]
        new_counts[disc] = new_tuplist
    counts = new_counts
    if loglog:
        new_counts = {}
        for d,tuplist in counts.items():
            tuplist_loglog = []
            for tup in tuplist:
                tup0_loglog = log(max(1, tup[0]))
                tup1_loglog = log(max(1, tup[1]))
                tuplist_loglog.append((tup0_loglog, tup1_loglog))
            new_counts[d] = tuplist_loglog
        counts = new_counts
    model = get_model(loglog=loglog)
    fits = {}
    plots = {}
    for d in disc_list:
        fits[d] = get_fit_from_counts(counts, d, model, spacing=spacing)
        if loglog:
            scale=None
        else:
            scale='loglog'
        plots[d] = make_plot(counts, d, model, fits[d], show_plot=False, scale=scale)
    if show_plots:
        for d in disc_list:
            print d,':',fits[d]
        show(sum(flatten(plots.values())))
    to_ret = (plots, fits)
    return to_ret


def make_counts_by_deg(min_level=10000, max_level=2000000, degs=None, disc_list=None, starting_counts=None, level_data=None, loglog=False):
    if degs is None:
        degs = list(range(1,7))
    if starting_counts is None:
        starting_counts = {1:329, 2:212, 3:76, 4:28, 5:20, 6:11}
        #starting_counts = {1:329, 2:212, 3:76, 4:10, 5:20, 6:11} # This was for checking how the exponent in the fit depends on the starting count for deg 4
        #starting_counts = {1:329, 2:210, 3:70, 4:23, 5:2, 6:0} # This was with only discriminants that appear > 10k
    if level_data is None:
        level_data = read_level_to_deg_sign_disc()
    counts_by_deg = {d:[(previous_prime(min_level), starting_counts[d])] for d in degs}
    for level in primes(min_level, max_level):
        level_counts = {d:0 for d in degs}
        if level in level_data:
            for d,sgn,disc in level_data[level]:
                if (d in level_counts) and ((disc_list is None) or (disc in disc_list)):
                    level_counts[d] += 1
        for d in degs:
            prev_count = counts_by_deg[d][-1][1]
            counts_by_deg[d].append((level, prev_count + level_counts[d]))
    if loglog:
        new_counts_by_deg = {}
        for d,tuplist in counts_by_deg.items():
            tuplist_loglog = []
            for tup in tuplist:
                tup0_loglog = log(max(1, tup[0]))
                tup1_loglog = log(max(1, tup[1]))
                tuplist_loglog.append((tup0_loglog, tup1_loglog))
            new_counts_by_deg[d] = tuplist_loglog
        counts_by_deg = new_counts_by_deg
    return counts_by_deg


def get_model(ab_vars=None, x_var=None, loglog=False):
    if ab_vars is None:
        ab_vars = var('a,b')
    a_var, b_var = ab_vars
    if x_var is None:
        x_var = var('x')
    if loglog:
        model = log(a_var*log_integral(exp(x_var)**b_var))
    else:
        model = a_var*log_integral(x_var**b_var)
    return model


def get_fit_from_counts(counts_by_deg, deg, model, spacing=1, ab_vars=None, x_var=None):
    if ab_vars is None:
        ab_vars = var('a,b')
    if x_var is None:
        x_var = var('x')
    data = counts_by_deg[deg][::spacing]
    fit = find_fit(data, model, solution_dict=True, parameters=list(ab_vars), variables=[x_var])#, initial_guess=(6.0,1.0))
    return fit


def make_plot(counts_by_deg, deg, model, fit, min_level=None, max_level=None, show_plot=False, scale='loglog'):
    x_var = [v for v in model.variables() if v not in fit][0]
    if min_level is None:
        min_level = min([tmp[0] for tmp in counts_by_deg[deg]])
    if max_level is None:
        max_level = max([tmp[0] for tmp in counts_by_deg[deg]])
    pt_plot = points(counts_by_deg[deg], scale=scale)
    fit_plot = plot(model.subs(fit), (x_var, min_level, max_level), scale=scale, color='red')
    if show_plot:
        show(pt_plot + fit_plot)
    return (pt_plot, fit_plot)


def make_deg_count_plots(fits=None, loglog=True, spacing=100, show_plot=True, dpi=800):
    degs = [1,2,3,4]
    counts_by_deg = make_counts_by_deg(degs=degs, loglog=loglog)
    counts_by_deg_nologlog = make_counts_by_deg(degs=degs, loglog=False)
    model = get_model(loglog=loglog)
    model_nologlog = get_model(loglog=False)
    
    a,b,x = var('a,b,x')
    if fits is None:
        if (spacing == 100) and loglog:
            fits = {}
            fits[1] = {b: 0.8321696787655243, a: 0.9687240229795021}
            fits[2] = {b: 0.6201631982508222, a: 3.389847515714529}
            fits[3] = {b: 0.33070917752443674, a: 7.364485074837851}
            fits[4] = {b: 0.10722004335315158, a: 11.961534915746284}
        elif (spacing == 30) and loglog:
            # With smaller spacing these stay the same to the precision displayed in the figure
            fits = {}
            fits[1] = {b: 0.8321628521354324, a: 0.9688019264853368}
            fits[2] = {b: 0.6200648636693188, a: 3.393770915430526}
            fits[3] = {b: 0.3307994171148687, a: 7.357964360072556}
            fits[4] = {b: 0.10717010807768236, a: 11.96841930729677}
        else:
            fits = {}
            for d in degs:
                fits[d] = get_fit_from_counts(counts_by_deg, d, model, spacing=spacing)
    
    colour_list = [hue(0.0, s=0.8), hue(0.12, s=0.8), hue(0.26, s=0.8), hue(0.48, s=0.8)]
    size = 12
    pts1 = points(counts_by_deg_nologlog[1], scale='loglog', color=colour_list[0], dpi=dpi, size=size, legend_label='Degree 1', legend_color=colour_list[0])
    pts2 = points(counts_by_deg_nologlog[2], scale='loglog', color=colour_list[1], dpi=dpi, size=size, legend_label='Degree 2', legend_color=colour_list[1])
    pts3 = points(counts_by_deg_nologlog[3], scale='loglog', color=colour_list[2], dpi=dpi, size=size, legend_label='Degree 3', legend_color=colour_list[2])
    pts4 = points(counts_by_deg_nologlog[4], scale='loglog', color=colour_list[3], dpi=dpi, size=size, legend_label='Degree 4', legend_color=colour_list[3])
    pts1.set_legend_options(markerscale=1.8)
    pts2.set_legend_options(markerscale=1.8)
    pts3.set_legend_options(markerscale=1.8)
    pts4.set_legend_options(markerscale=1.8)
    
    fp1 = plot(model_nologlog.subs(fits[1]), (x, 10**4, 2*10**6), ymin=1, scale='loglog', color='black', thickness=2)
    fp2 = plot(model_nologlog.subs(fits[2]), (x, 10**4, 2*10**6), ymin=1, scale='loglog', color='black', thickness=2)
    fp3 = plot(model_nologlog.subs(fits[3]), (x, 10**4, 2*10**6), ymin=1, scale='loglog', color='black', thickness=2)
    fp4 = plot(model_nologlog.subs(fits[4]), (x, 10**4, 2*10**6), ymin=1, scale='loglog', color='black', thickness=2)
    
    t1str = '$' + str(fits[1][a]+0.005)[:4] + '\,\mathrm{li}(X^{' + str(fits[1][b]+0.0005)[:5] + '})$'
    t2str = '$' + str(fits[2][a]+0.005)[:4] + '\,\mathrm{li}(X^{' + str(fits[2][b]+0.0005)[:5] + '})$'
    t3str = '$' + str(fits[3][a]+0.005)[:4] + '\,\mathrm{li}(X^{' + str(fits[3][b]+0.0005)[:5] + '})$'
    t4str = '$' + str(fits[4][a]+0.05)[:4] + '\,\mathrm{li}(X^{' + str(fits[4][b]+0.0005)[:5] + '})$'
    
    positions = [(0.5*10**5, 2.23*10**3), (2.1*10**5, 0.63*10**3), (0.9*10**5, 1.9*10**2), (2.3*10**5, 22)]
    t1 = text(t1str, positions[0], rotation=0.0, color='black', fontsize=14)
    t2 = text(t2str, positions[1], rotation=0.0, color='black', fontsize=14)
    t3 = text(t3str, positions[2], rotation=0.0, color='black', fontsize=14)
    t4 = text(t4str, positions[3], rotation=0.0, color='black', fontsize=14)

    plt = pts1+pts2+pts3+pts4+fp1+fp2+fp3+fp4+t1+t2+t3+t4
    if show_plot:
        plt.show()
    
    return plt



###############################
### Poisson log-likelihoods ###
###############################

'''
sage: al_25 = [0.5+tmp*0.005 for tmp in range(1401)]
sage: bl_25 = [-0.28 - tmp*0.0001 for tmp in range(2001)]
sage: rat_25 = 1.2
sage: true_count_25 = 1900+986
sage: true_count_25
2886
sage: sum(3*p**(-0.36) for p in primes(10**4, 2*10**6))
3750.12496876978
sage: 37/29.0
1.27586206896552
sage: rat_25 = 1.3
sage: sum(2*p**(-0.36) for p in primes(10**4, 2*10**6))
2500.08331251318
sage: lpthresh_25 = {}
sage: time tmplp = dp.get_poisson_logprobs_with_countcheck(al_25, bl_25, true_count_25, rat_25, deg_list=[2], disc_list=[5], toupdate=lpthresh_25)
CPU times: user 26min 41s, sys: 273 ms, total: 26min 41s
Wall time: 26min 43s
sage: time diyptsthresh_25 = dp.logprob_dict_to_diy_contour_plot_list(lpthresh_25, dpi=400); show(sum(diyptsthresh_25))
Launched png viewer for Graphics object consisting of 6 graphics primitives
CPU times: user 17.5 s, sys: 79.6 ms, total: 17.6 s
Wall time: 17.7 s
sage: save(lpthresh_25, 'poisson_logprobs_deg2_disc5')
sage: al_1 = [0.4 + 0.001*tmp for tmp in range(1501)]
sage: bl_1 = [-0.11 - 0.00005*tmp for tmp in range(2001)]
sage: sum(0.7*p**(-0.15) for p in primes(10**4, 2*10**6))
13882.7481590128
sage: true_count_1 = 8843 + 6406
sage: true_count_1
15249
sage: 13882.8/15249
0.910407239819005
sage: 1/0.91
1.09890109890110
sage: rat_1 = 1.2
sage: lpthresh_1 = {}
sage: time tmplp = dp.get_poisson_logprobs_with_countcheck(al_1, bl_1, true_count_1, rat_1, deg_list=[1], disc_list=None, toupdate=lpthresh_1)
CPU times: user 25min 51s, sys: 175 ms, total: 25min 51s
Wall time: 25min 52s
sage: time diyptsthresh_1 = dp.logprob_dict_to_diy_contour_plot_list(lpthresh_1, dpi=400); show(sum(diyptsthresh_1))
Launched png viewer for Graphics object consisting of 6 graphics primitives
CPU times: user 12.9 s, sys: 100 ms, total: 13 s
Wall time: 12.9 s
sage: save(lpthresh_1, 'poisson_logprobs_deg1')
sage: al_28 = [0.01*tmp for tmp in range(1,2000)]
sage: bl_28 = [-0.3 - tmp*0.0005 for tmp in range(1001)]
sage: al_28 = [0.01*tmp for tmp in range(1,2001)]
sage: true_count_28 = 342
sage: sum(p**(-0.5) for p in primes(10**4, 2*10**6))
205.642571218496
sage: 205*2.3
471.500000000000
sage: 205*2.3/342
1.37865497076023
sage: rat_28 = 1.7
sage: lpthresh_28 = {}
sage: time tmplp = dp.get_poisson_logprobs_with_countcheck(al_28, bl_28, true_count_28, rat_28, deg_list=[2], disc_list=[8], toupdate=lpthresh_28)
CPU times: user 14min 15s, sys: 520 ms, total: 14min 16s
Wall time: 14min 18s
sage: time diyptsthresh_28 = dp.logprob_dict_to_diy_contour_plot_list(lpthresh_28, dpi=400); show(sum(diyptsthresh_28))
Launched png viewer for Graphics object consisting of 6 graphics primitives
CPU times: user 12 s, sys: 68 ms, total: 12.1 s
Wall time: 12.1 s
sage: save(lpthresh_28, 'poisson_logprobs_deg2_disc8')

sage: al_349 = [tmp*0.01 for tmp in range(1,2501)]
sage: bl_349 = [-0.4 - tmp*0.0002 for tmp in range(3001)]
sage: sum(p**(-0.6) for p in primes(10**4, 2*10**6))
57.9235970926669
sage: rat_349 = 2.0
sage: lpthresh_349 = {}
sage: true_count_349 = 120
sage: time tmplp = dp.get_poisson_logprobs_with_countcheck(al_349, bl_349, true_count_349, rat_349, deg_list=[3], disc_list=[49], toupdate=lpthresh_349)
CPU times: user 40min 26s, sys: 621 ms, total: 40min 26s
Wall time: 40min 28s

sage: time diyptsthresh_349 = dp.logprob_dict_to_diy_contour_plot_list(lpthresh_349, dpi=400); show(sum(diyptsthresh_349))
Launched png viewer for Graphics object consisting of 6 graphics primitives
CPU times: user 39.4 s, sys: 428 ms, total: 39.9 s
Wall time: 39.9 s
sage: save(lpthresh_349, 'poisson_logprobs_deg3_disc49')

sage: bl_349 = [-0.2 - tmp*0.0002 for tmp in range(1001)]
sage: time tmplp = dp.get_poisson_logprobs_with_countcheck(al_349, bl_349, true_count_349, rat_349, deg_list=[3], disc_list=[49], toupdate=lpthresh_349)
CPU times: user 12min 2s, sys: 51.6 ms, total: 12min 2s
Wall time: 12min 2s

'''

def get_poisson_logprobs(a_list, b_list, deg_list=None, sign_list=None, disc_list=None, ldsd=None, toupdate=None, verbose=False):
    if ldsd is None:
        ldsd = read_level_to_deg_sign_disc()
    if toupdate is None:
        toupdate = {}
    
    level_list = list(primes(10**4, 2*10**6))
    counts = {}
    fac_term = 0
    for level in level_list:
        f_list = ldsd[level]
        count = 0
        for deg,sgn,disc in f_list:
            ok = True
            if (deg_list is not None) and (deg not in deg_list):
                ok = False
            if (sign_list is not None) and (sign not in sign_list):
                ok = False
            if (disc_list is not None) and (disc not in disc_list):
                ok = False
            if ok:
                count += 1
        if count > 0:
            if count not in counts:
                counts[count] = []
            counts[count].append(level)
            fac_term -= log(float(factorial(count)))
    if verbose:
        for c, ll in sorted(counts.items()):
            if len(ll) < 10:
                llcstr = str(ll)
            else:
                llcstr = str(ll[:5]) + ' ... ' + str(ll[-5:])
            print c, llcstr
    
    logN_vals = {level:log(float(level)) for level in level_list}
    logN_sums = {}
    for c in counts:
        c_term = 0
        for level in counts[c]:
            c_term += logN_vals[level]
        logN_sums[c] = c_term

    for b in b_list:
        level_sum = 0
        for level in level_list:
            level_sum -= float(level)**b
        
        for a in a_list:
            loga = log(float(a))
            
            count_term = 0
            for c in counts:
                c_term = logN_sums[c]
                c_term *= b
                c_term += loga * len(counts[c])
                c_term *= c
                count_term += c_term
            
            level_term = a*level_sum

            logprob = count_term + level_term + fac_term
            toupdate[(a,b)] = logprob
            if verbose:
                print 'count_term:', count_term
                print 'level_term:', level_term
                print 'fac_term:', fac_term
    
    return toupdate


def get_poisson_logprobs_with_countcheck(a_list, b_list, true_total, thresh, deg_list=None, sign_list=None, disc_list=None, ldsd=None, toupdate=None, verbose=False):
    if ldsd is None:
        ldsd = read_level_to_deg_sign_disc()
    if toupdate is None:
        toupdate = {}
    
    level_list = list(primes(10**4, 2*10**6))
    counts = {}
    fac_term = 0
    for level in level_list:
        f_list = ldsd[level]
        count = 0
        for deg,sgn,disc in f_list:
            ok = True
            if (deg_list is not None) and (deg not in deg_list):
                ok = False
            if (sign_list is not None) and (sign not in sign_list):
                ok = False
            if (disc_list is not None) and (disc not in disc_list):
                ok = False
            if ok:
                count += 1
        if count > 0:
            if count not in counts:
                counts[count] = []
            counts[count].append(level)
            fac_term -= log(float(factorial(count)))
    
    logN_vals = {level:log(float(level)) for level in level_list}
    logN_sums = {}
    for c in counts:
        c_term = 0
        for level in counts[c]:
            c_term += logN_vals[level]
        logN_sums[c] = c_term

    for b in b_list:
        level_sum = 0
        for level in level_list:
            level_sum -= float(level)**b
        
        for a in a_list:
            level_term = a*level_sum
            ratio = -level_term/true_total
            if (ratio < thresh) and (1/ratio < thresh):
                loga = log(float(a))
                count_term = 0
                for c in counts:
                    c_term = logN_sums[c]
                    c_term *= b
                    c_term += loga * len(counts[c])
                    c_term *= c
                    count_term += c_term
                logprob = count_term + level_term + fac_term
                toupdate[(a,b)] = logprob
    
    return toupdate


def logprob_dict_to_diy_contour_plot_list(logprobs, maxtup=None, levels=None, color_list=None, a_range=None, b_range=None, max_pts_len=None, **kwargs):
    if maxtup is None:
        maxk, maxval = sorted(logprobs.items(), key = lambda tup: -tup[1])[0]
    if levels is None:
        num_levels = 6
        #levels = [maxval - (i+1)*log(10.0) for i in range(num_levels)]
        levels = [maxval - i*log(10.0) - 0.000001 for i in range(num_levels)]
    levels = sorted(levels, reverse=True)
    if color_list is None:
        phi = 0.5*(1 + sqrt(5.0))
        color_list = [hue(i/phi) for i in range(1+len(levels))]
    point_plots = []
    for i,lvl in enumerate(levels):
        goodpts = []
        for k,v in logprobs.items():
            good = True
            if a_range is not None:
                if (k[0] < a_range[0]) or (a_range[1] < k[0]):
                    good = False
            if b_range is not None:
                if (k[1] < b_range[0]) or (b_range[1] < k[1]):
                    good = False
            if v < lvl:
                good = False
            if (i>0) and (v > levels[i-1]):
                good = False
            if good:
                goodpts.append(k)
        if max_pts_len is not None:
            random.shuffle(goodpts)
            goodpts = goodpts[:max_pts_len]
        if i == 0:
            kwargs_copy = dict(kwargs)
            kwargs_copy['size'] = 10
            pts = points(goodpts, color=color_list[i], faceted=True, markeredgecolor='black', **kwargs_copy)
        else:
            if 'size' not in kwargs:
                pts = points(goodpts, color=color_list[i], size=2, **kwargs)
            else:
                pts = points(goodpts, color=color_list[i], **kwargs)
        point_plots.append(pts)
    point_plots.reverse()
    return point_plots



########################
### Poisson p-values ###
########################


def get_collision_p_values(verbose=True):
    actual_collision_counts = get_collision_counts()
    expected_collision_counts = get_expected_collision_counts()
    p_values = {}
    for d in expected_collision_counts:
        p_values_d = {}
        for k in expected_collision_counts[d]:
            actual = actual_collision_counts[d].get(k,0)
            expected = expected_collision_counts[d][k]
            tail_mass = get_tail_estimate(actual, expected)
            p_values_d[k] = tail_mass
            if verbose:
                if tail_mass == 0:
                    print 'tail_mass == 0 at d,k =', d, k
                    print 'actual:', actual
                    print 'expected:', expected
                    log_chernoff_bound = actual * (1 + log(expected)) - expected - actual * log(actual) # from wikipedia
                    print 'chernoff_bound:', RealField(400)(exp(log_chernoff_bound))
                    if actual > expected:
                        print 'using sf:'
                        print poisson.sf(actual-0.5,expected)
        p_values[d] = p_values_d
    return p_values


def le_cam(disc, k_list=None, verbose=False):
    if disc not in [1,5,8,49]:
        raise NotImplementedError
    if k_list is None:
        k_list = list(range(13))
    a,b = get_MLE_abs()[disc]
    
    means = {}
    square_sums = {}
    for k in k_list:
        if verbose:
            print k
        prob_sum = 0
        prob_squared_sum = 0
        for p in primes(10**4, 2*10**6):
            mu = a*p**b
            prob = poisson.pmf(k,mu) # from scipy
            prob_sum += prob
            prob_squared_sum += prob**2
        means[k] = prob_sum
        square_sums[k] = prob_squared_sum
    
    le_cam_errors = {}
    for k in square_sums:
        err = 2*square_sums[k]
        mu = means[k]
        if mu > 1:
            err /= mu
        le_cam_errors[k] = err

    to_ret = (means, square_sums, le_cam_errors)
    return to_ret


def get_MLE_abs():
    MLE_abs = {}
    MLE_abs[1] = (0.897, -0.1616)
    MLE_abs[5] = (2.795, -0.3747)
    MLE_abs[8] = (1.41, -0.4870)
    MLE_abs[49] = (0.74, -0.5186)
    return MLE_abs


def get_collision_counts():
    count_tups = [(1, {1: 11922, 2: 1038, 3: 302, 4: 60, 5: 15, 6: 1, 7: 2, 10: 1}), (5, {1: 2772, 2: 54, 3: 2}), (8, {1: 336, 2: 3}), (12, {1: 17}), (13, {1: 46}), (17, {1: 1}), (21, {1: 4}), (49, {1: 120}), (81, {1: 13}), (148, {1: 6}), (169, {1: 9}), (229, {1: 19, 2: 1}), (257, {1: 7}), (321, {1: 1}), (725, {1: 6}), (1957, {1: 2}), (2777, {1: 2}), (8768, {1: 1}), (14641, {1: 1}), (70601, {1: 1}), (371293, {1: 1})]
    collision_counts = {d:counts for (d,counts) in count_tups}
    for d in collision_counts:
        zero_count = 147704 - sum(collision_counts[d].values())
        collision_counts[d][0] = zero_count
    
    return collision_counts


def get_tail_estimate(count, mu):
    if count < mu:
        tail_mass = poisson.cdf(count+0.5, mu) # prob Pois(mu) <= count
    else:
        tail_mass = 1 - poisson.cdf(count-0.5, mu) # prob Pois(mu) >= count
    return tail_mass


def get_expected_collision_counts():
    '''
    sage: time lc1 = dp.le_cam(1, k_list=list(range(13)))
    CPU times: user 4min 2s, sys: 292 ms, total: 4min 2s
    Wall time: 4min 4s
    sage: pprint.pprint(lc1)
    ({0: 133240.09894865425,
    1: 13708.878323305627,
    2: 727.4793946774062,
    3: 26.75089255641614,
    4: 0.7731501681826703,
    5: 0.018874946339583806,
    6: 0.0004075298913455458,
    7: 8.014618889541047e-06,
    8: 1.4612845405373093e-07,
    9: 2.493607869859711e-09,
    10: 4.0003547990300574e-11,
    11: 6.045090514718649e-13,
    12: 8.614373766881187e-15},
    {0: 120232.7038023499,
    1: 1303.3656343976145,
    2: 4.107706777178623,
    3: 0.007106356439964928,
    4: 8.777336272387832e-06,
    5: 8.540475963391319e-09,
    6: 6.684456872266939e-12,
    7: 4.222546687387003e-15,
    8: 2.1679255527014676e-18,
    9: 9.150030464207973e-22,
    10: 3.216210295851556e-25,
    11: 9.535714564222549e-29,
    12: 2.4131118820178782e-32},
    {0: 1.8047525444826198,
    1: 0.19014912871198836,
    2: 0.011292984536009151,
    3: 0.0005312986417165722,
    4: 1.7554672544775664e-05,
    5: 1.7080951926782638e-08,
    6: 1.3368913744533879e-11,
    7: 8.445093374774006e-15,
    8: 4.335851105402935e-18,
    9: 1.8300060928415946e-21,
    10: 6.432420591703112e-25,
    11: 1.9071429128445099e-28,
    12: 4.8262237640357564e-32})
    
    sage: time lc5 = dp.le_cam(5, k_list=list(range(6)))
    CPU times: user 1min 52s, sys: 212 ms, total: 1min 52s
    Wall time: 1min 53s
    sage: pprint.pprint(lc5)
    ({0: 144852.82698006407,
    1: 2816.524929730339,
    2: 34.266365209678874,
    3: 0.3774514048628623,
    4: 0.0042259051042054345,
    5: 4.7180516896349534e-05},
    {0: 142070.19514335907,
    1: 66.31780408737566,
    2: 0.023984149346120923,
    3: 9.366071717329156e-06,
    4: 2.914975260261294e-09,
    5: 6.695187853939492e-13},
    {0: 1.9615798753159581,
    1: 0.047091934736558566,
    2: 0.0013998653898281783,
    3: 1.8732143434658312e-05,
    4: 5.829950520522588e-09,
    5: 1.3390375707878983e-12})
    
    sage: time lc8 = dp.le_cam(8, k_list=list(range(5)))
    CPU times: user 1min 32s, sys: 156 ms, total: 1min 32s
    Wall time: 1min 33s
    sage: pprint.pprint(lc8)
    ({0: 147362.25194286028,
    1: 341.1398976105717,
    2: 0.606937206286318,
    3: 0.0012196164704385767,
    4: 2.6993327007471023e-06},
    {0: 147021.71776553502,
    1: 1.2065889890006667,
    2: 1.6021409137221777e-05,
    3: 2.31992692135361e-10,
    4: 2.4627640924055026e-15},
    {0: 1.9953782712623407,
    1: 0.007073866161372004,
    2: 3.2042818274443554e-05,
    3: 4.63985384270722e-10,
    4: 4.925528184811005e-15})
    
    sage: time lc49 = dp.le_cam(49, k_list=list(range(5)))
    CPU times: user 1min 32s, sys: 228 ms, total: 1min 32s
    Wall time: 1min 34s
    sage: pprint.pprint(lc49)
    ({0: 147584.03180780742,
    1: 119.88777837506241,
    2: 0.08034920593291892,
    3: 6.455161789255781e-05,
    4: 5.7474425887644546e-08},
    {0: 147464.22431414505,
    1: 0.1603117908596192,
    2: 3.433603677923705e-07,
    3: 7.847923403286851e-13,
    4: 1.2946726049438657e-18},
    {0: 1.9983764165784765,
    1: 0.0026743641934558575,
    2: 6.86720735584741e-07,
    3: 1.5695846806573703e-12,
    4: 2.5893452098877315e-18})
'''
    expected_collision_counts = {}
    expected_collision_counts[1] = {0: 133240.09894865425,
                                    1: 13708.878323305627,
                                    2: 727.4793946774062,
                                    3: 26.75089255641614,
                                    4: 0.7731501681826703,
                                    5: 0.018874946339583806,
                                    6: 0.0004075298913455458,
                                    7: 8.014618889541047e-06,
                                    8: 1.4612845405373093e-07,
                                    9: 2.493607869859711e-09,
                                    10: 4.0003547990300574e-11,
                                    11: 6.045090514718649e-13,
                                    12: 8.614373766881187e-15}
    expected_collision_counts[5] = {0: 144852.82698006407,
                                    1: 2816.524929730339,
                                    2: 34.266365209678874,
                                    3: 0.3774514048628623,
                                    4: 0.0042259051042054345,
                                    5: 4.7180516896349534e-05}
    expected_collision_counts[8] = {0: 147362.25194286028,
                                    1: 341.1398976105717,
                                    2: 0.606937206286318,
                                    3: 0.0012196164704385767,
                                    4: 2.6993327007471023e-06}
    expected_collision_counts[49] = {0: 147584.03180780742,
                                     1: 119.88777837506241,
                                     2: 0.08034920593291892,
                                     3: 6.455161789255781e-05,
                                     4: 5.7474425887644546e-08}
    return expected_collision_counts


def get_random_poisson_logprobs(true_a, true_b, hypothesis_a, hypothesis_b, num_trials, toupdate=None, use_real_data=False, real_data_info=None):
    # This was used to check that the log-likelihoods were similar to log-likelihoods of random data generated by the model.
    toupdate_was_None = toupdate is None
    if toupdate is None:
        toupdate = []
    if use_real_data:
        ldsd = read_level_to_deg_sign_disc()
        deg_list, sign_list, disc_list = real_data_info
    
    true_means = {}
    level_sum = 0
    log_facs = {}
    for N in primes(10**4, 2*10**6):
        true_means[N] = true_a*float(N)**true_b
        level_sum -= hypothesis_a*float(N)**hypothesis_b
        log_facs[N] = log(hypothesis_a*float(N)**hypothesis_b)

    if use_real_data:
        num_trials = 1
    for _ in range(num_trials):
        logprob = level_sum
        for N in primes(10**4, 2*10**6):
            if use_real_data:
                c = 0
                for (tmp_deg, tmp_sign, tmp_disc) in ldsd[N]:
                    good = True
                    if not ((deg_list is None) or (tmp_deg in deg_list)):
                        good = False
                    if not ((sign_list is None) or (tmp_sign in sign_list)):
                        good = False
                    if not ((disc_list is None) or (tmp_disc in disc_list)):
                        good = False
                    if good:
                        c += 1
            else:
                c = numpy.random.poisson(true_means[N])
            logprob += c*log_facs[N]
            if c > 1:
                logprob -= log(float(factorial(c)))
        toupdate.append(logprob)
    
    if toupdate_was_None:
        to_ret = toupdate
    else:
        to_ret = None
    return to_ret




################
### AL signs ###
################


def get_random_walk_prob(level_to_deg_sign_disc, deg, weight_power_gen, weight_power_eval_list, disc_list=None, omit_SN=False, min_level=10**4, level_list=None, use_real_data=False, normalize=False):
    # This, with use_real_data=True, was used to give likelihoods of each beta generating the data.
    if weight_power_gen == 0:
        def expected_fn_gen(tmp1,tmp2):
            return 0.5
    else:
        expected_fn_gen = get_expected_fn_power(weight_power_gen)
    
    expected_fn_eval_dict = {}
    for wt in weight_power_eval_list:
        expected_fn_eval_dict[wt] = get_expected_fn_power(wt)
    
    prob_dict = {wt:RR(1) for wt in weight_power_eval_list}
    #SN_levels = {}
    SN_levels = set(get_setzer_neumann_levels())
    if level_list is None:
        level_list = sorted(level_to_deg_sign_disc.keys())
    for level in level_list:
        if level > min_level:
            v = level_to_deg_sign_disc[level]
            sn_ok = (not omit_SN) or (level not in SN_levels) #(not (sqrt(level - 64) in ZZ))
            for (tmp_deg, tmp_sign, tmp_disc) in v:
                if sn_ok or (not (tmp_deg, tmp_sign) == (1,-1)):
                    if (tmp_deg == deg) and ((disc_list is None) or (tmp_disc in disc_list)):
                        plus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == -1])
                        
                        if use_real_data:
                            random_sign = tmp_sign
                        else:
                            plus_prob = 1 - expected_fn_gen(plus_dim, minus_dim) # expected_fn is prob of minus
                            r = random.random()
                            if r < plus_prob:
                                random_sign = 1
                            else:
                                random_sign = -1
                        
                        for wt in prob_dict:
                            plus_prob = 1 - expected_fn_eval_dict[wt](plus_dim, minus_dim)
                            if random_sign == 1:
                                prob_dict[wt] *= plus_prob
                            else:
                                prob_dict[wt] *= 1-plus_prob
                            if normalize:
                                prob_dict[wt] *= 2
                else:
                    sn_ok = True
                    #print 'Removing SN. Level = '+str(level)+', v = '+str(v)
                    #if level not in SN_levels:
                    #    SN_levels[level] = []
                    #SN_levels[level].append((tmp_deg, tmp_sign, tmp_disc))
    to_ret = prob_dict #prob #(prob, SN_levels)
    return to_ret


def get_random_walk_prob_SW(weight_power_eval_list, min_index=0, max_index=99, omit_SN=False, min_level=10**4, max_level=None, normalize=False):
    expected_fn_eval_dict = {}
    for wt in weight_power_eval_list:
        expected_fn_eval_dict[wt] = get_expected_fn_power(wt)
    
    prob_dict = {wt:RR(1) for wt in weight_power_eval_list}
    SN_levels = set()
    for index in range(min_index, max_index+1):
        print index
        for C in SteinWatkinsPrimeData(index):
            level = C.conductor
            if (level > min_level) and ((max_level is None) or (level < max_level)):
                sgn = (C.rank % 2) * 2 - 1
                sn_ok = (not omit_SN) or (sgn == 1) or (level in SN_levels) or (not (sqrt(level - 64) in ZZ))
                if not sn_ok:
                    SN_levels.add(level)
                else:
                    plus_dim = level/24.0 - 0.5*sqrt(float(level))
                    minus_dim = level/24.0 + 0.5*sqrt(float(level))
                    for wt in prob_dict:
                        plus_prob = 1 - expected_fn_eval_dict[wt](plus_dim, minus_dim)
                        if sgn == 1:
                            prob_dict[wt] *= plus_prob
                        else:
                            prob_dict[wt] *= 1-plus_prob
                        if normalize:
                            prob_dict[wt] *= 2
    to_ret = prob_dict #prob #(prob, SN_levels)
    return to_ret



def get_random_walk(level_to_deg_sign_disc, deg, weight_power, disc_list=None, omit_SN=False, min_level=10**4, level_list=None, shift_by_EV=False, use_real_data=False):
    if weight_power == 0:
        def expected_fn(tmp1,tmp2):
            return 0.5
    else:
        expected_fn = get_expected_fn_power(weight_power)
    count_diffs = {}
    diff_tally = 0
    if level_list is None:
        level_list = sorted(level_to_deg_sign_disc.keys())
    for level in level_list:
        if level > min_level:
            v = level_to_deg_sign_disc[level]
            sn_ok = (not omit_SN) or (not (sqrt(level - 64) in ZZ))
            for (tmp_deg, tmp_sign, tmp_disc) in v:
                if sn_ok or (not (tmp_deg, tmp_sign) == (1,-1)):
                    if (tmp_deg == deg) and ((disc_list is None) or (tmp_disc in disc_list)):
                        plus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == 1])
                        minus_dim = sum([tmp[0] for tmp in level_to_deg_sign_disc[level] if tmp[1] == -1])
                        plus_prob = 1 - expected_fn(plus_dim, minus_dim) # expected_fn is prob of minus
                        r = random.random()
                        if r < plus_prob:
                            random_sign = 1
                        else:
                            random_sign = -1
                        if use_real_data:
                            random_sign = tmp_sign
                        diff_tally += random_sign
                        if shift_by_EV:
                            diff_tally += 1 - 2*plus_prob # -( 1*(1-expected_fn) + (-1)*expected_fn )
                        count_diffs[level] = diff_tally
                else:
                    sn_ok = True
    return count_diffs


def p_value_pts(rwprobs_d, deltasum, alpha_bounds=None, beta_bounds=None, spacing=1):
    # Prob(Prob(random seq from beta | alpha) > Prob(data | alpha))
    # deltasum_bydeg = {1: 0.01992139147761528, 2: 0.007743806493261199, 3: 0.0006077076080877266, 4: 0.00020439259658961724, 5: 5.237898425222656e-06, 6: 1.4814473468199136e-06}
    # deltasum_bydeg.update({-1: 0.01972085834725267, -2: 0.007743088761474358, -3: 0.0006076809880740136, -4: 0.00020439259658961724, -5: 5.237898425222656e-06, -6: 1.4814473468199136e-06})
    prob_vals = {}
    for alpha in rwprobs_d.keys()[::spacing]:
        if (alpha_bounds is None) or ((alpha_bounds[0] < alpha) and (alpha < alpha_bounds[1])):
            for beta in rwprobs_d.keys()[::spacing]:
                if (beta_bounds is None) or ((beta_bounds[0] < beta) and (beta < beta_bounds[1])):
                    mu = (144*alpha*beta - 72*alpha**2) * deltasum
                    var = (144*alpha**2) * deltasum
                    zscore = (log(rwprobs_d[alpha]) - mu)/sqrt(2*var + 10**(-15)) # This is zscore/sqrt2
                    pval = 0.5*(1 + erf(zscore))
                    if pval > 0.5:
                        pval = 1 - pval
                    prob_vals[(alpha,beta)] = pval
    return prob_vals


def get_p_values(rw_probs=None, dataprobs=None, deltasum=None, deltasum_bydeg=None):
    # This compares expressions for mean and variance + central limit theorem to random walks generated by the model.
    if rw_probs is None:
        rw_probs = load('random_walk_probs_beta0.sobj')
        # rw_probs = {1:{1:probs1_0_1, -0.17:probs1_0_m}, 2:{1:probs2_0_1, -0.38:probs2_0_m}, 3:{1:probs3_0_1, -0.67:probs3_0_m}, -1:{1:probs1nosn_0_1, -0.17:probs1nosn_0_m}}
    bd_vals = {1:-0.17, 2:-0.38, 3:-0.67}
    dataprobs_was_none = dataprobs is None
    if dataprobs is None:
        ldsd = read_level_to_deg_sign_disc()
        dataprobs = {d:{tmp:get_random_walk_prob(ldsd, abs(d), 0, [tmp], use_real_data=True, normalize=True, omit_SN=(d==-1))[tmp] for tmp in [0,1.0,bd_vals[abs(d)]]} for d in [1,-1,2,3]}
    #ldsd = None
    if deltasum is None:
        #ldsd = read_level_to_deg_sign_disc()
        #deltasum = 0.0
        #for p in primes(10**4, 2*10**6):
        #    plus_dim = sum([tmp[0] for tmp in ldsd[p] if tmp[1] == 1])
        #    minus_dim = sum([tmp[0] for tmp in ldsd[p] if tmp[1] == -1])
        #    deltasum += (plus_dim - minus_dim)**2 / float(p)**2
        deltasum = 0.16172544967832755
    if deltasum_bydeg is None:
        deltasum_bydeg = {1: 0.01992139147761528, 2: 0.007743806493261199, 3: 0.0006077076080877266, 4: 0.00020439259658961724, 5: 5.237898425222656e-06, 6: 1.4814473468199136e-06}
        deltasum_bydeg.update({-1: 0.01972085834725267, -2: 0.007743088761474358, -3: 0.0006076809880740136, -4: 0.00020439259658961724, -5: 5.237898425222656e-06, -6: 1.4814473468199136e-06})
    
    zscores = {d:{} for d in [1,-1,2,3]}
    mustd_dict = {d:{} for d in [1,-1,2,3]}
    for d in [1,-1,2,3]:
        print '\n'
        print 'Degree:', d
        # actually generated by beta
        # what is prob it was generated by alpha?
        for alpha in [0, 1, bd_vals[abs(d)]]:
            for beta in [0, 1, bd_vals[abs(d)]]:
                print ''
                print 'alpha, beta: ' + '%.2f' % alpha + ', ' + '%.2f' % beta
                mu = 144*alpha*beta - 72*alpha**2
                var = 144*alpha**2
                ds = deltasum_bydeg[d]
                mu *= ds
                var *= ds
                mustd_dict[d][(alpha,beta)] = (mu,var**0.5)
                print 'mean, std: ' + '%.4f' % mu + ' +- ' + '%.4f' % (var**0.5)
                dataprb_alpha = dataprobs[d][alpha]
                rho = dataprb_alpha
                if var != 0:
                    zscore = (log(rho) - mu) / var**0.5
                else:
                    zscore = None
                zscores[d][(alpha,beta)] = zscore
                print 'rho, log(rho): ' + '%.4f' % rho + ', ' + '%.4f' % log(rho)
                print 'zscore: ' + str(zscore)[:4]
                if zscore is not None:
                    cdf_val = 0.5 * (1 + erf(zscore/sqrt(2.0)))
                    print 'CDF value: ' + '%.4f' % cdf_val
                if beta == 0 and alpha != 0:
                    rwpl = rw_probs[d][alpha]
                    print 'Number of samples:', len(rwpl)
                    mu_empirical = sum([log(tmp) for tmp in rwpl])/len(rwpl)
                    std_empirical = std([log(tmp) for tmp in rwpl])
                    print 'Empirical: ' + '%.4f' % mu_empirical + ' +- ' + '%.4f' % std_empirical
                    success_ratio = len([tmp for tmp in rwpl if tmp < rho])/float(len(rwpl))
                    print 'Success ratio: ' + '%.4f' % success_ratio
    to_ret = [zscores, mustd_dict]
    if dataprobs_was_none:
        to_ret.append(dataprobs)
    to_ret = tuple(to_ret)
    return to_ret



####################
### Lang-Trotter ###
####################


def build_form_dicts(deg, num_levels_at_once=40, min_level=10000, max_level=None, print_info=True, prime_only=None, run_remove_shorter_dupes=True, num_eigs=None, verbose=True, write_filename=None, dicts_to_update=None):
    al_data = load('al_data.sobj')
    levels = []
    for k,v in sorted(al_data.items()):
        for d,s in v:
            if (d == deg) and ((min_level is None) or (k > min_level)) and ((max_level is None) or (k < max_level)):
                levels.append(k)
    dirlist1 = get_dirs_external_hdd()
    dirlist2 = get_dirs_external_hdd_1M_2M()
    if dicts_to_update is None:
        formlens_dict = {}
        pia_counts_dict = {}
        pia_counts_dict_all = {}
        rational_ans_count_dict = {}
        rational_ans_count_dict_all = {}
    else:
        formlens_dict = dicts_to_update[0]
        pia_counts_dict = dicts_to_update[1]
        pia_counts_dict_all = dicts_to_update[2]
        rational_ans_count_dict = dicts_to_update[3]
        rational_ans_count_dict_all = dicts_to_update[4]
    num_chunks = len(levels)/num_levels_at_once + 1
    for i in range(num_chunks):
        level_chunk = levels[i*num_levels_at_once : (i+1)*num_levels_at_once]
        hev_dict = read_eigs_from_dirs(dirlist1 + dirlist2, num_eigs=num_eigs, levels=level_chunk, degs=[deg], check_level_bounds_from_filename=True, verbose=verbose)
        if run_remove_shorter_dupes:
            new_hev_dict = {}
            for k in hev_dict:
                if type(hev_dict[k]) == type([]):
                    f_list = hev_dict[k]
                    new_f_list = remove_shorter_dupes(f_list, verbose=verbose)
                    new_hev_dict[k] = new_f_list
                else:
                    for cond in hev_dict[k]:
                        f_list = hev_dict[k][cond]
                        new_f_list = remove_shorter_dupes(f_list, verbose=verbose)
                        if k not in new_hev_dict:
                            new_hev_dict[k] = {}
                        new_hev_dict[k][cond] = new_f_list
            hev_dict = new_hev_dict
        if prime_only is None:
            if print_info:
                print 'pia_counts:'
            if write_filename is not None:
                wf = open(write_filename,'a')
                wf.write('pia_counts:\n')
                wf.close()
            chunk_pia_counts_dict, chunk_formlens_dict = hev_dict_to_pia_counts(hev_dict, prime_only=True, print_info=print_info, with_formlens=True, write_filename=write_filename)
            if print_info:
                print 'pia_counts_all:'
            if write_filename is not None:
                wf = open(write_filename,'a')
                wf.write('pia_counts_all:\n')
                wf.close()
            chunk_pia_counts_dict_all = hev_dict_to_pia_counts(hev_dict, prime_only=False, print_info=print_info, with_formlens=False, write_filename=write_filename)
            if print_info:
                print 'rational_ans:'
            if write_filename is not None:
                wf = open(write_filename,'a')
                wf.write('rational_ans:\n')
                wf.close()
            chunk_rational_ans_count_dict = hev_dict_to_rational_ans(hev_dict, prime_only=True, print_info=print_info, with_formlens=False, write_filename=write_filename)
            if print_info:
                print 'rational_ans_all:'
            if write_filename is not None:
                wf = open(write_filename,'a')
                wf.write('rational_ans_all:\n')
                wf.close()
            chunk_rational_ans_count_dict_all = hev_dict_to_rational_ans(hev_dict, prime_only=False, print_info=print_info, with_formlens=False, write_filename=write_filename)
        else:
            if prime_only:
                if print_info:
                    print 'pia_counts:'
                if write_filename is not None:
                    wf = open(write_filename,'a')
                    wf.write('pia_counts:\n')
                    wf.close()
                chunk_pia_counts_dict, chunk_formlens_dict = hev_dict_to_pia_counts(hev_dict, prime_only=True, print_info=print_info, with_formlens=True, write_filename=write_filename)
                chunk_pia_counts_dict_all = {}
                if print_info:
                    print 'rational_ans:'
                if write_filename is not None:
                    wf = open(write_filename,'a')
                    wf.write('rational_ans:\n')
                    wf.close()
                chunk_rational_ans_count_dict = hev_dict_to_rational_ans(hev_dict, prime_only=True, print_info=print_info, with_formlens=False, write_filename=write_filename)
                chunk_rational_ans_count_dict_all = {}
            else:
                chunk_pia_counts_dict = {}
                if print_info:
                    print 'pia_counts_all:'
                if write_filename is not None:
                    wf = open(write_filename,'a')
                    wf.write('pia_counts_all:\n')
                    wf.close()
                chunk_pia_counts_dict_all, chunk_formlens_dict = hev_dict_to_pia_counts(hev_dict, prime_only=False, print_info=print_info, with_formlens=True, write_filename=write_filename)
                chunk_rational_ans_count_dict = {}
                if print_info:
                    print 'rational_ans_all:'
                if write_filename is not None:
                    wf = open(write_filename,'a')
                    wf.write('rational_ans_all:\n')
                    wf.close()
                chunk_rational_ans_count_dict_all = hev_dict_to_rational_ans(hev_dict, prime_only=False, print_info=print_info, with_formlens=False, write_filename=write_filename)
        formlens_dict.update(chunk_formlens_dict)
        pia_counts_dict.update(chunk_pia_counts_dict)
        pia_counts_dict_all.update(chunk_pia_counts_dict_all)
        rational_ans_count_dict.update(chunk_rational_ans_count_dict)
        rational_ans_count_dict_all.update(chunk_rational_ans_count_dict_all)
    if dicts_to_update is None:
        to_ret = (formlens_dict, pia_counts_dict, pia_counts_dict_all, rational_ans_count_dict, rational_ans_count_dict_all)
    else:
        to_ret = None
    return to_ret
            
            
def hev_dict_to_rational_ans(hev_dict, prime_only=True, print_info=False, with_formlens=False, write_filename=None):
    rational_ans_dict = {}
    formlens_dict = {}
    for level in sorted(hev_dict):
        if print_info:
            print ''
            print ''
            print level
        if write_filename is not None:
            wf = open(write_filename, 'a')
            wf.write('\n\n'+str(level)+'\n')
            wf.close()
        rational_ans_dict[level] = []
        formlens_dict[level] = []
        if type(hev_dict[level]) == type([]):
            vals = {None: hev_dict[level]}
        else:
            vals = hev_dict[level]
        for cond in vals:
            if print_info:
                print ''
                print 'AL sign: ' + str(cond if (cond is None) else -cond)
            if write_filename is not None:
                wf = open(write_filename, 'a')
                wf.write('AL sign: '+str(cond if (cond is None) else -cond)+'\n')
                wf.close()
            for f in vals[cond]:
                if print_info:
                    print 'len(f):', len(f)
                nf = get_K(f)
                if print_info:
                    print 'Number field:'
                    print '\t', nf
                    print '\t', nf.galois_group('pari')
                if write_filename is not None:
                    wf = open(write_filename, 'a')
                    wf.write('len(f): '+str(len(f))+'\n')
                    wf.write('Number field:\n')
                    wf.write('\t' + str(nf)+'\n')
                    wf.write('\t' + str(nf.galois_group('pari'))+'\n')
                    wf.close()
                formlens_dict[level].append((len(f), nf))
                rational_ans = get_rational_ans(f, prime_only=prime_only)
                rational_ans_dict[level].append(rational_ans)
                if print_info or (write_filename is not None):
                    if write_filename is not None:
                        wf = open(write_filename,'a')
                    if len(rational_ans.keys()) < 50:
                        for n,an in sorted(rational_ans.items()):
                            if print_info:
                                print an, '\t', n#,'\t', ZZ(n).factor() if n != 0 else 0
                            if write_filename is not None:
                                wf.write(str(an)+'\t'+str(n)+'\n')
                    else:
                        counts = {}
                        for n,an in rational_ans.items():
                            if an not in counts:
                                counts[an] = []
                            counts[an].append(n)
                        already_printed = set()
                        for n,an in sorted(rational_ans.items()):
                            if an not in already_printed:
                                if len(counts[an]) == 1:
                                    if print_info:
                                        print an, '\t', n#,'\t', ZZ(n).factor() if n != 0 else 0
                                    if write_filename is not None:
                                        wf.write(str(an)+'\t'+str(n)+'\n')
                                elif len(counts[an]) < 10:
                                    if print_info:
                                        print an,'\t', sorted(counts[an])
                                    if write_filename is not None:
                                        wf.write(str(an)+'\t'+str(sorted(counts[an]))+'\n')
                                else:
                                    if print_info:
                                        print an,'\t', sorted(counts[an])[:25],'...',sorted(counts[an])[-25:]
                                    if write_filename is not None:
                                        wf.write(str(an)+'\t'+str(sorted(counts[an])[:25])+'...'+str(sorted(counts[an])[-25:])+'\n')
                                already_printed.add(an)
                    if write_filename is not None:
                        wf.close()
    if with_formlens:
        to_ret = (rational_ans_dict, formlens_dict)
    else:
        to_ret = rational_ans_dict
    return to_ret


def hev_dict_to_pia_counts(hev_dict, prime_only=True, print_info=False, with_formlens=False, write_filename=None):
    pia_counts_dict = {}
    formlens_dict = {}
    for level in sorted(hev_dict):
        if print_info:
            print ''
            print ''
            print level
        if write_filename is not None:
            wf = open(write_filename, 'a')
            wf.write('\n\n'+str(level)+'\n')
            wf.close()
        pia_counts_dict[level] = []
        formlens_dict[level] = []
        if type(hev_dict[level]) == type([]):
            vals = {None: hev_dict[level]}
        else:
            vals = hev_dict[level]
        for cond in vals:
            if print_info:
                print ''
                print 'AL sign: ' + str(cond if (cond is None) else -cond)
            if write_filename is not None:
                wf = open(write_filename, 'a')
                wf.write('AL sign: '+str(cond if (cond is None) else -cond)+'\n')
                wf.close()
            for f in vals[cond]:
                if print_info:
                    print 'len(f):', len(f)
                nf = get_K(f)
                if print_info:
                    print 'Number field:'
                    print '\t', nf
                    print '\t', nf.galois_group('pari')
                if write_filename is not None:
                    wf = open(write_filename, 'a')
                    wf.write('len(f): '+str(len(f))+'\n')
                    wf.write('Number field:\n')
                    wf.write('\t' + str(nf)+'\n')
                    wf.write('\t' + str(nf.galois_group('pari'))+'\n')
                    wf.close()
                formlens_dict[level].append((len(f), nf))
                an_to_ns = get_an_to_ns(f, prime_only=prime_only)
                counts = {}
                for an in an_to_ns:
                    c = len(an_to_ns[an])
                    if c not in counts:
                        counts[c] = [an]
                    elif counts[c] in ZZ:
                        counts[c] += 1
                    else:
                        counts[c].append(an)
                        if len(counts[c]) > 10:
                            counts[c] = len(counts[c])
                pia_counts_dict[level].append({tmp1:tmp2 if tmp2 in ZZ else len(tmp2) for tmp1,tmp2 in counts.items()})
                if print_info or (write_filename is not None):
                    if write_filename is not None:
                        wf = open(write_filename,'a')
                    for c, an_list in sorted(counts.items()):
                        if an_list in ZZ:
                            if print_info:
                                print c,'\t', an_list
                            if write_filename is not None:
                                wf.write(str(c)+'\t'+str(an_list)+'\n')
                        else:
                            if print_info:
                                print c,'\t', len(an_list)
                            if write_filename is not None:
                                wf.write(str(c)+'\t'+str(len(an_list))+'\n')
                            for an in an_list:
                                nl = sorted(an_to_ns[an])
                                #nl = [ZZ(tmp).factor() for tmp in nl]
                                if len(nl) < 50:
                                    nl_str = str(nl)
                                else:
                                    nl_str = str(nl[:25]) + ' ... ' + str(nl[-25:])
                                if print_info:
                                    print '\t\t' + str(an) + '\t' + nl_str
                                if write_filename is not None:
                                    wf.write('\t\t' + str(an) + '\t' + nl_str + '\n')
                    if write_filename is not None:
                        wf.close()
    if with_formlens:
        to_ret = (pia_counts_dict, formlens_dict)
    else:
        to_ret = pia_counts_dict
    return to_ret


def get_an_to_ns(f, prime_only=False, counts_only=False):
    an_to_ns = {}
    if type(f) == type([]):
        if prime_only:
            n_gen = primes(len(f))
        else:
            n_gen = range(len(f))
    else:
        n_gen = f.keys()
        if prime_only:
            n_gen = [tmp for tmp in n_gen if ZZ(tmp).is_prime()]
    for n in n_gen:
        an = f[n]
        if not counts_only:
            if an not in an_to_ns:
                an_to_ns[an] = []
            an_to_ns[an].append(n)
        else:
            if an not in an_to_ns:
                an_to_ns[an] = 0
            an_to_ns[an] += 1
    return an_to_ns


def get_subfield_hits(f, deg_check=None):
    K = get_K(f)
    if not ((deg_check is None) or (K.degree() == deg_check)):
        return None
    S = [tmp[0] for tmp in K.subfields()[:-1]]
    subfield_hits = {}
    for n,an in enumerate(f):
        if n not in [0,1]:
            for L in S:
                conj_list = an.galois_conjugates(L)
                if len(conj_list) > 0:
                    if n not in subfield_hits:
                        subfield_hits[n] = [an, [(L, conj_list)]]
                    else:
                        subfield_hits[n][1].append((L, conj_list))
    return subfield_hits


def get_rational_ans(f, prime_only=False):
    rational_ans = {}
    if prime_only:
        n_gen = primes(len(f))
    else:
        n_gen = range(len(f))
    for n in n_gen:
        an = f[n]
        if an in ZZ:
            rational_ans[n] = an
    return rational_ans


def make_pia_counts_by_D(form_dicts, countvals=None, normalize=False):
    # output of build_form_dicts
    # (formlens_galoisgroups, pia_counts, pia_counts_all, ratans, ratans_all)
    if countvals is None:
        countvals = list(range(2,10))
    max_c = max(countvals)
    pia_counts_by_D = {c:{} for c in countvals}
    for level in sorted(form_dicts[0].keys()):
        for i in range(len(form_dicts[0][level])):
            formlen = form_dicts[0][level][i][0]
            K = form_dicts[0][level][i][1]
            if hasattr(K, 'number_field'):
                K = K.number_field()
            D = K.disc()
            for c in pia_counts_by_D:
                v = form_dicts[1][level][i].get(c,0)
                if normalize:
                    if K.degree() == 2:
                        v /= float(prime_pi(formlen))
                    elif K.degree() == 3:
                        v /= float(prime_pi(formlen))**(0.5)
                if D not in pia_counts_by_D[c]:
                    pia_counts_by_D[c][D] = []
                pia_counts_by_D[c][D].append(v)
    return pia_counts_by_D


def form_dicts_to_pia_points(form_dicts, countvals=None, by_disc=True, normalize=False):
    if countvals is None:
        countvals = list(range(2,10))
    counts = {} #{c:[] for c in countvals}
    for level in sorted(form_dicts[1].keys()):
        for i in range(len(form_dicts[1][level])):
            formlen = form_dicts[0][level][i][0]
            K = form_dicts[0][level][i][1]
            if hasattr(K, 'number_field'):
                K = K.number_field()
            D = K.disc()
            if not by_disc:
                D = None
            for c in countvals:
                v = form_dicts[1][level][i].get(c,0)
                if normalize:
                    if K.degree() == 2:
                        v /= float(prime_pi(formlen))
                    elif K.degree() == 3:
                        v /= float(prime_pi(formlen))**(0.5)
                if D not in counts:
                    counts[D] = {c:[] for c in countvals}
                counts[D][c].append((level,v))
    if not by_disc:
        counts = counts[None]
    return counts

  
def form_dicts_to_ratans_points(form_dicts, by_disc=True, normalize=False):
    counts = {} #{c:[] for c in countvals}
    for level in sorted(form_dicts[3].keys()):
        for i in range(len(form_dicts[3][level])):
            formlen = form_dicts[0][level][i][0]
            K = form_dicts[0][level][i][1]
            if hasattr(K, 'number_field'):
                K = K.number_field()
            D = K.disc()
            if not by_disc:
                D = None
            v = len(form_dicts[3][level][i])
            if normalize:
                if K.degree() == 2:
                    v /= float(prime_pi(formlen))**0.5
                elif K.degree() == 3:
                    v /= log(float(prime_pi(formlen)))
            if D not in counts:
                counts[D] = []
            counts[D].append((level,v))
    if not by_disc:
        counts = counts[None]
    return counts


def ratans_points_to_plot(counts, colour_list=None, size=15, show_plot=True):
    plt_list = []
    by_disc = type(counts) == type({})
    if not by_disc:
        counts = {None: counts}
    if colour_list is None:
        if len(counts) == 6:
            colour_list = [hue(0.0), hue(0.55), hue(0.28), hue(0.67), hue(0.12), hue(0.78)]
        elif len(counts) == 7:
            colour_list = [hue(0.58), hue(0.70), hue(0.8), hue(0.0), hue(0.12), hue(0.26), hue(0.46)]
            #colour_list = [hue(0.65), hue(0.26), hue(0.79), hue(0.13), hue(0.95), hue(0.55), hue(0.07)]
        else:
            colour_list = [hue(i/float(len(counts))) for i in range(len(counts))]
    discs = sorted(counts.keys())
    for i,D in enumerate(discs):
        colour = colour_list[i]
        faceted = len(counts[D]) < 10
        if faceted:
            markeredgecolor = 'black'
        else:
            markeredgecolor = None
        label_str = 'D = '+str(D)
        if by_disc:
            pts = points([(tmp1,tmp2 + 0.0 * (random.random() - 0.5)) for tmp1,tmp2 in counts[D]], color=colour, size=size, legend_color=colour, legend_label=label_str, faceted=faceted, markeredgecolor=markeredgecolor)
        else:
            pts = points(counts[D], color=colour, size=size)
        plt_list.append(pts)
    if show_plot:
        show(sum(plt_list), figsize=12)
    return plt_list


def pia_counts_by_D_to_histo(pia_counts_by_D, colour_list=None, bins=None, show_histo=True):
    # despite name, expects dict {D:count_list}, not {c:{D:count_list}}
    discs = sorted(pia_counts_by_D.keys())
    datalists = [pia_counts_by_D[D] for D in discs]
    max_len = max([len(dl) for dl in datalists])
    datalists_padded = [dl + [-1]*(max_len - len(dl)) for dl in datalists]
    if colour_list is None:
        if len(discs) == 6:
            colour_list = [hue(0.0), hue(0.55), hue(0.28), hue(0.67), hue(0.12), hue(0.78)]
        elif len(discs) == 7:
            colour_list = [hue(0.59), hue(0.72), hue(0.8), hue(0.0), hue(0.15), hue(0.3), hue(0.52)]
        else:
            colour_list = [hue(i/float(len(discs))) for i in range(len(discs))]
    if bins is None:
        max_val = max([max(dl) for dl in datalists])
        if max_val < 30:
            bins = list(range(max_val+2))
        else:
            bins = [max_val * i/20.0 for i in range(21)]
    h = histogram(datalists_padded, bins=bins, stacked=True, color=colour_list, label=['D = '+str(D) for D in discs])
    if show_histo:
        h.show(show_legend=True, legend_handlelength=1.0)
    return h


def make_c22_plot(counts, marker_list=None, colour_list=None, size=15, dpi=400, show_plot=True):
    if marker_list is None:
        marker_list = ['o', 's', '+', 'd', 'x', '*', 'p', 'D', 'H', 'h']
    plt_list = []
    if not by_disc:
        counts = {None: counts}
    if colour_list is None:
        if len(counts) == 6:
            colour_list = [hue(0.0), hue(0.55), hue(0.28), hue(0.67), hue(0.12), hue(0.78)]
        elif len(counts) == 7:
            colour_list = [hue(0.58), hue(0.70), hue(0.8), hue(0.0), hue(0.12), hue(0.26), hue(0.46)]
        else:
            colour_list = [hue(i/float(len(counts))) for i in range(len(counts))]
    discs = sorted(counts.keys())
    for i,D in enumerate(discs):
        colour = colour_list[i]
        faceted = len(counts[D]) < 20
        if faceted:
            markeredgecolor = 'black'
        else:
            markeredgecolor = None
        label_str = 'D = '+str(D)
        pts = points([(tmp1,tmp2 + 0.0 * (random.random() - 0.5)) for tmp1,tmp2 in counts[D]], color=colour, size=size, legend_color=colour, legend_label=label_str, faceted=faceted, markeredgecolor=markeredgecolor, dpi=dpi)
        plt_list.append(pts)
    if show_plot:
        show(sum(plt_list), figsize=12)
    return plt_list



#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    if '-profile' in sys.argv:
        cProfile.run('main()')
    else:
        main()
