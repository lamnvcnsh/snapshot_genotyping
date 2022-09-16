from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from scipy.signal import find_peaks,peak_prominences,peak_widths


class FSA:
    def __init__(self) -> None:
        pass


def loadFSA(aFSA):
    """Load FSA file

    Args:
        aFSA (str): A path to a FSA file.

    Returns:
        str: Name of FSA file.
        dict: A dictionary data loaded from FSA file.
    """

    return(os.path.basename(aFSA), SeqIO.read(aFSA, 'abi'))



def generateReferenceRange(reference=[15, 20, 25, 35, 50, 62, 80, 110, 120]):
    """Generate all range from reference data

    Args:
        reference (list, required): A list contain all bp locations of refernces.
                                    Defaults to [15, 20, 25, 35, 50, 62, 80, 110, 120].

    Returns:
        list:   The list contant all generated range of the references
                Example: [(15, 20), (20,25)]
    """

    referenceRange = []
    
    for i,_ in enumerate(reference):
        if i < len(reference)-1:
            referenceRange.append((reference[i], reference[i+1]))

    return referenceRange






print(help(generateReferenceRange))



#####################
def extract_raw_intensity(data):
    
    channels = ['DATA1', 'DATA2', 'DATA3', 'DATA4', 'DATA105']

    # DATA105 reference channel (LIZ120)

    intensity = {}

    for channel in channels:
        intensity[channel] = data[channel]
    
    return intensity



def generate_range(lists):
    range = []
    for i,_ in enumerate(lists):
        # print(i)
        if i < len(lists)-1:
            range.append((lists[i], lists[i+1]))

    return range



def liz120_info(ref_size =[15, 20, 25, 35, 50, 62, 80, 110, 120]):

    return ref_size, generate_range(ref_size)



def finding_peaks(intensity, min_width=1, min_height=500, min_prominence = 50):
    
    peaks, heights = find_peaks(intensity, height=min_height)
    heights = heights['peak_heights']

    if min_prominence:
        prominence = peak_prominences(intensity, peaks)
        conditions =  np.where(prominence[0] > min_prominence)
        peaks = peaks[conditions]
        heights = heights[conditions]
    
    if min_width:
        width = peak_widths(intensity, peaks)
        conditions = np.where(width[0] > min_width)
        peaks = peaks[conditions]
        heights = heights[conditions]

    return peaks, heights


def union(*list):
    final_list = []
    for l in list:
        final_list = final_list + l
        # print(l)
    return final_list


def generate_base_size_from_ref(ref_intensity, range=(0,200)):
    
    size, ranges = liz120_info()
    peaks, peak_heights = finding_peaks(ref_intensity, min_height=800)
    total_points = np.max(peaks) - np.min(peaks)
    size_range = np.max(size) - np.min(size)

    # print(peaks, total_points, size_range, size)
    
    # average points is total point divide to size range
    avarage_points = total_points/size_range

    # simulate first base range and last base range
    # since we need to extend from ref size to both size includes begin and end
    # default range from 0, 200 bp

    first_range = np.min(size) - range[0]
    
    first_range_points = int(round(avarage_points * first_range))
    first_range_bases = np.linspace(range[0], np.min(size), first_range_points).tolist()

    last_range = range[1] - np.max(size) 
    last_range_points = int(round(avarage_points * last_range))
    last_range_bases = np.linspace(np.max(size), range[1], last_range_points).tolist()

    # print(first_range, last_range, first_range_points, last_range_points)
    # simulate the reference base
    peak_ranges = generate_range(peaks)

    # expected that number of range from references and ranges form detected peaks are identical
    if len(ranges) != len(peak_ranges):
        raise ValueError('Reference range and detected peak range are difference!')

    base_ref_range = []

    for ref, peak in zip(ranges, peak_ranges):

        singe_range = np.linspace(ref[0], ref[1], peak[1] - peak[0]+ 1)
        base_ref_range = base_ref_range + singe_range.tolist()
    
    # gathering all base range 
    full_base = first_range_bases + base_ref_range + last_range_bases
    full_base = list(np.round(np.array(full_base), 2))
    data_points_index = (np.min(peaks) - first_range_points + 1, np.max(peaks) + last_range_points)

    # create data

    data = {}
    data['reference_size'] = size
    data['full_base'] = np.sort(np.unique(np.array(full_base)))
    data['intensity'] = np.array(ref_intensity[data_points_index[0]:data_points_index[1]])
    data['point_index'] = data_points_index
    data['detected_peaks'] = peaks
    data['peak_heights'] = peak_heights


    return data
###########################

print(generateReferenceRange([1,2,3]))

