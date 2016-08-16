import numpy as np
import matplotlib.pyplot as plt
from sunpy.io import fits
import sunpy.map
import os
import datetime
import time
from sunpy.lightcurve import LightCurve
from scipy.optimize import curve_fit
from sunpy.physics.solar_rotation import mapcube_solar_derotate
import scipy
import batman
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.interpolate import interp1d

#Fixes a bug in sunpy.
#This is fixed by https://github.com/sunpy/sunpy/pull/1853
sunpy.time.time.TIME_FORMAT_LIST.append("%Y.%m.%d_%H:%M:%S_TAI")

personal_data_directory = '/home/gtaylor/'
aia_data_directory = '/data/SDO/AIA/level1/'
hmi_data_directory = '/data/SDO/HMI/intensitygram/'

#The angle to rotate the picture of the sun so that Venus goes horizontal through the image
transit_rotate = {"venus":7.4,"mercury":2.0}

no_rotate = 0
#Splits are a dictionary of two tuples of two tuples each that are each in the format (w,[x,y,z,...]). 
#That is, it breaks the image into w slices and then uses slices x,y,z,...
#The first tuple is for AIA images, the second tuple is for HMI images
#The section of the sun that contains the transit.
#(64,[13,14])
inside_split = {"venus":((32,[8]), (16,[3])),"mercury":((256,[159,160,161]), (256,[166,167,168]))}
#The section of the sun that doesn't contain the transit
outside_split = {"venus":(((32,range(8) + range(9,32)), (16,range(3) + range(4,16)))),"mercury":((256,range(159) + range(162,256)), ((256,range(166) + range(169,256))))}
no_split = {"venus":(((1,[0]), (1,[0]))),"mercury":(((1,[0]), (1,[0])))}

venus_spike_start_time = 9000 #in seconds, start of spike in 1700Å
venus_spike_end_time = 15000 #in seconds, end of spike in 1700Å


ingress_start_time = {"venus":4020,"mercury":5010} #in seconds, start of ingress
ingress_end_time = {"venus":5150,"mercury":5280} #in seconds, end of ingress

egress_start_time = {"venus":26100,"mercury":31440} #in seconds, start of egress
egress_end_time = {"venus":27350,"mercury":31680} #in seconds, end of egress

transit_start_time = ingress_end_time
transit_end_time = egress_start_time

def timeSinceStartOfTransit(planet, time):
    '''
    `planet` is the planet transit to use

    Given a `time` since the start of the recorded period, 
        returns the difference between the given time 
        and when the transit actually started
    '''
    return time - transit_start_time[planet]

def timeThroughTransit(planet, time):
    '''
    `planet` is the planet transit to use

    Given a `time` since the start of the recorded period, 
        returns    the start time if the transit has yet to start,
                   the end time   if the transit has already finished,
                or the current time
    '''
    if time < transit_start_time[planet]:
        return transit_start_time[planet]
    if time < transit_end_time[planet]:
        return time
    return transit_end_time[planet]        

def checkTimeBlock(planet, time_block, time):
    '''
    `planet` is the planet transit to use

    `time_block` is an array of strings to limit the list to (see section "Time")
   
    `time` is the time to check

    Checks if `time` is in the `time_block`
    '''
    if time_block == []:
        return True
    
    doIAppend = False
    
    if "transit" in time_block:
        if time > transit_start_time[planet] and time < transit_end_time[planet]:
            doIAppend = True
            
    if "not transit" in time_block:
        if time < ingress_start_time[planet] or time > egress_end_time[planet]:
            doIAppend = True

    if "ingress" in time_block:
        if time > ingress_start_time[planet] and time < ingress_end_time[planet]:
            doIAppend = True
    if "egress" in time_block:
        if time > egress_start_time[planet] and time < egress_end_time[planet]:
            doIAppend = True

    if "not venus spike" in time_block:
        if time > venus_spike_start_time and time < venus_spike_end_time:
            doIAppend = False
    
    return doIAppend

def getDataFromTimeBlock(planet, time_block, data, times):
    '''
    `planet` is the planet transit to use

    `time_block` is an array of strings to limit the list to (see section "Time")

    `data` is an array of generic data points that correspond with the time points in `times`

    `data` and `times` must be of same length    

    Returns the data points that match the `time_block`
    '''
    good_data = []
    for i, time in enumerate(times):
        if checkTimeBlock(planet, time_block, time):
            good_data.append(data[i])
            
    return good_data

#data about Venus

venus_radius = 3760.4 #miles
sun_radius = 432168.6 #miles
venus_semi_major_axis = 67237909. #miles

#apparent venus radius in stellar radii
#((tangent of angular radius of venus on june 5th 2012)*(1 astronomical unit - radius of sun - radius of earth) / (radius of sun))
venus_apparant_radius = 0.023 #stellar radii

venus_longitude_perihelion = 131.53298 #degrees

venus_orbital_period_days = 224.701 #days
venus_orbital_period = venus_orbital_period_days * 24 * 60 * 60 #seconds

venus_eccentricity = 0.0067
venus_orbital_inclination = 3.39 #degrees

venus_time_inferior_conjunc = (transit_start_time["venus"] + transit_end_time["venus"] ) / 2 #seconds


'''
Functions that are of the form nth * x^n + n-1th * x^(n-1) + … + one * x + zero

They work for `func_one` and `func_two` for graphLightCurve
'''

def cubic(x,three,two,one,zero):
    return np.multiply(three,np.power(x,3)) + np.multiply(two,np.square(x)) + np.multiply(one,x) + zero

def quadratic(x,two,one,zero):
    return np.multiply(two,np.square(x)) + np.multiply(one,x) + zero

def linear(x,one,zero):
    return np.multiply(one,x) + zero

def quartic(x,four,three,two,one,zero):
    return np.multiply(four,np.power(x,4)) + np.multiply(three,np.power(x,3)) + np.multiply(two,np.square(x)) + np.multiply(one,x) + zero

def quintic(x,five,four,three,two,one,zero):
    return np.multiply(five,np.power(x,5)) + np.multiply(four,np.power(x,4)) + np.multiply(three,np.power(x,3)) + np.multiply(two,np.square(x)) + np.multiply(one,x) + zero


wavelengths = ["0094","0131","0171","0193","0211","0304","0335","1600","1700","cont"]
#"cont" = HMI intensity continuum, all other are AIA wavelengths

def checkIfInWavelengths(wavelength):
    '''
    Raises an exception if the given `wavelength` 
        is not in the list of wavelengths
    '''
    if wavelength not in wavelengths:
        raise Exception('Wavelength ' + wavelength + 
                        ' not an accepted wavelength. Try one of ' 
                        + str(wavelengths))
        
def getWavelengthFromFile(filepath):
    '''
    `filepath` is a path to a FITS file
    '''
    if "intensity" in filepath:
        return "cont"
    else:
        return filepath[-9:-5]


def getPlanetFromFile(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns either "venus" or "mercury"
    '''
    if "/2012/" in filepath:
        return "venus"
    elif "/2016/" in filepath:
        return "mercury"
    else:
        raise Exception("Can't parse filepath " + filepath + " to get a planet")


def readData(planet, wavelength):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    Returns the saved data for that wavelength in list form, sorted by time
    '''
    checkIfInWavelengths(wavelength)
    
    data = sorted(readFromFile(personal_data_directory + planet + 'data' + wavelength + '.txt').items(), key=lambda tup: tup[1][1])
        
    good_data = []
    
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, pixels_out_of_split) in data:
        if np.isfinite(count):
            good_data.append((filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                                         pixels_out_of_split)))
    
    return good_data

def writeToFile(filepath, dictionary):
    '''
    `filepath` is the filepath to write to

    `dictionary` is the dictionary to save

    Returns the dictionary
    '''
    with open(filepath, 'a') as target:
        target.truncate(0)
        target.write(str(dictionary))
        
    return dictionary
        
def readFromFile(filepath):
    '''
    `filepath` is the filepath to read from
    
    Returns the file, evaluated. 
    
    This will turn, for example, a dictionary in the file 
        into a dictionary in python
    '''
    while not os.path.exists(filepath):
        print("Waiting for file " + filepath + " to exist")
        time.sleep(1)
    
    try:
        with open(filepath, 'r') as target:
            s = target.read()
            return eval(s)
    except:
        return {}


def showFancyImage(filepath, subplot=None, new_figure=True):
    '''
    `filepath` is a path to a FITS file

    `subplot` is the argument passed into matplotlib's subplot function. If `subplot == None` then subplot is not used.

    If `new_figure` is True, it creates a new matplotlib figure for the image

    It makes a sunpy map of the given file, and puts the map on screen
    
    Returns the map
    '''

    if new_figure:
        plt.figure()

    if subplot != None:
        plt.subplot(subplot)

    smap = sunpy.map.Map(filepath)
    smap.plot()
    return smap

def showArray(array, subplot=None, new_figure=True):
    '''
    Graphs the given `array` in greyscale 

    `subplot` is the argument passed into matplotlib's subplot function. If `subplot == None` then subplot is not used.

    If `new_figure` is True, it creates a new matplotlib figure for the graphed array

    Returns the array
    '''
    
    if new_figure:
        plt.figure()
        
    if subplot != None:
        plt.subplot(subplot)
        
    plt.matshow(array, cmap='Greys')
    
    return array
    

def showHistogramOfPixelIntensities(array, boxes=50, subplot=None, new_figure=True):
    '''
    Shows a histogram of elements of the given `array`

    `boxes` is the number of separate bars to sort the elements of the array into
    
    `subplot` is the argument passed into matplotlib's subplot function. If `subplot == None` then subplot is not used.

    If `new_figure` is True, it creates a new matplotlib figure for the histogram

    Note: uses lots of ram
        
    Returns the array, with all nan elements set to zero
    '''
            
    if new_figure:
        plt.figure()
        
    if subplot != None:
        plt.subplot(subplot)

    array = array[~np.isnan(array)]
    
    x = np.random.random_integers(0,10000000)
    plt.figure(x)
    plt.hist(array, boxes)
    plt.show()
    
    return array

def showImage(filepath, split=no_split, rotate=no_rotate, replace_with_zero=False, subplot=None, new_figure=True):
    '''
    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    If `replace_with_zero` is true, it will replace the removed elements of the array with zeros, otherwise it graphs a truncated array.

    `subplot` is the argument passed into matplotlib's subplot function. If `subplot == None` then subplot is not used.

    If `new_figure` is True, it creates a new matplotlib figure for the graphed array

    Graphs the array from the filepath, applying the split and rotate
    '''
            
    if new_figure:
        plt.figure()
        
    if subplot != None:
        plt.subplot(subplot)

    array = getRotatedAndSplitArrayFromFitsFile(filepath, split, rotate, replace_with_zero=replace_with_zero)
    
    showArray(array)
    
    return array


def getArrayFromFitsFile(filepath):
    '''
    `filepath` is a path to a FITS file

    If the file doesn't exist, wait for it to.

    Returns the array from the given file
    '''
    
    while not os.path.exists(filepath):
        print("Waiting for file " + filepath + " to exist")
        time.sleep(1)

    try:
        return np.array(fits.read(filepath)[1][0])
    except:
        print("Couldn't get array from " + filepath)
        time.sleep(1)
        return getArrayFromFitsFile(filepath)
    
def getHeaderFromFitsFile(filepath):
    '''
    `filepath` is a path to a FITS file

    If the file doesn't exist, wait for it to.

    Returns the header from the given file
    '''
    while not os.path.exists(filepath):
        print("Waiting for file " + filepath + " to exist")
        time.sleep(1)

    try:
        return fits.get_header(filepath)
    except:
        print("Couldn't get header from " + filepath)
        time.sleep(1)
        return getHeaderFromFitsFile(filepath)
        
def findFileByTimeSinceStartOfTransit(planet, time, wavelength):
    '''
    `planet` is the planet transit to use

    `time` is the time since the start of the transit to get the file from

    `wavelength` is the wavelength to use

    Returns the filepath to the FITS file closest to that time
    '''
    
    checkIfInWavelengths(wavelength)

    file_and_time_list = readListOfFilesAndTimeSinceStartOfTransit(planet, wavelength)
    filepath = None
    plusminus = 1
    
    if file_and_time_list == None:
        return None
    
    while filepath == None:
        filepath = next((i for i, v in enumerate(file_and_time_list) if v[1] < time + plusminus 
                     and v[1] > time - plusminus), None)
        if filepath == None:
            plusminus *= 2
    return file_and_time_list[filepath][0]

def getRotatedAndSplitArrayFromFitsFile(filepath, split=no_split, rotate=no_rotate, replace_with_zero=False):
    '''
    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    If `replace_with_zero` is true, it will replace the removed elements of the array with zeros, otherwise it returns a truncated array.

    Finds the array from the given file given file, 
        then applies the given split and rotate parameters as needed
        
    Returns the array created
    '''
    
    wavelength = getWavelengthFromFile(filepath)
    
    planet = getPlanetFromFile(filepath)
    
    if rotate == transit_rotate:
        rotate = transit_rotate[planet]
            
    split = split[planet]
        
    #Get the array
    full_array = getArrayFromFitsFile(filepath)
    if full_array == None:
        return []
    
    #Replace all the nans with zero
    full_array = np.nan_to_num(full_array)
    
    #Flips the array so that the orientation is consistent for all wavelengths
    
    
    if wavelength == "cont":
        full_array = np.fliplr(full_array)
        split = split[1]
    else:
        full_array = np.flipud(full_array)
        split = split[0]
            
    #Rotates the image, if neeeded
    if rotate != 0:
        full_array = scipy.ndimage.interpolation.rotate(full_array,rotate, reshape=False)
    
    #Splits the array into the different sections
    arrays = np.array_split(full_array, split[0])
    
    split_shape = arrays[0].shape
    
    #Goes through all the split arrays and isolates the wanted splits
    good_arrays = []    
    for i, array in enumerate(arrays):
        if i in split[1]:
            good_arrays.append(array)
        else:
            if replace_with_zero:
                good_arrays.append(np.zeros(split_shape))
         
    #Combines the wanted splits into a single image
    array = np.concatenate(good_arrays)
    return array

def getTotalData(filepath, split=no_split, rotate=no_rotate):
    '''
    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)
    
    It ignores the negative items in the array and items with a value 
        greater than the mean + std dev of the wavelength as a whole
    
    Returns the sum of the image data in the fits file
    '''
    
    wavelength = getWavelengthFromFile(filepath)
        
    planet = getPlanetFromFile(filepath)
    
    array = getRotatedAndSplitArrayFromFitsFile(filepath, split, rotate)
    
    #Remove the bottom of cut-off images
    array[np.where(array == -32768)] = 0
    #Cap pixels at the mean plus one standard deviation
    array[np.where(array > mean_stddev_skew_median_dict[planet][wavelength][0] 
                   + mean_stddev_skew_median_dict[planet][wavelength][1])] = (
        mean_stddev_skew_median_dict[planet][wavelength][0] 
        + mean_stddev_skew_median_dict[planet][wavelength][1])
    
    return np.nansum(array)

def getNumberOfPixels(filepath, split=no_split, rotate=no_rotate):
    '''
    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)
    
    Returns the number of non-nan pixels
    '''
    
    return np.count_nonzero(~np.isnan(getRotatedAndSplitArrayFromFitsFile(filepath, split, rotate)))

def getExpTime(filepath):
    '''
    `filepath` is a path to a FITS file

    Note: the HMI continuum doesn't have an exposure time, so returns a 1

    Returns the saved exposure time of the file. This may be inaccurate, depending on the processing already done on the file
    '''
    if getWavelengthFromFile(filepath) == "cont":
        return 1.0
    return getHeaderFromFitsFile(filepath)[1]["EXPTIME"]

def getActualExpTime(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the actual exposure time of the file, based on the time stamps in the file header
    '''
    return (sunpy.time.parse_time(getHeaderFromFitsFile(filepath)[1]["T_OBS"]) - 
            sunpy.time.parse_time(getHeaderFromFitsFile(filepath)[1]["DATE-OBS"])).microseconds / 1000. / 1000. * 2.0

def getCountsPerSecond(filepath, split=no_split, rotate=no_rotate):
    '''
    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    Returns the number of counts per second from the file, applying the split and rotation angle

    '''
    
    return getTotalData(filepath, split=split, rotate=rotate)/getExpTime(filepath)

def getTimeSinceStartOfTransit(planet, filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the time in seconds since the start of the recorded transit period
    '''
    
    planet = getPlanetFromFile(filepath)
    
    if planet == "venus":
        start_time = datetime.datetime(2012,6,5,21,0,0,0)
    else:
        start_time = datetime.datetime(2016,5,9,10,0,0,0)
        
    current_time = sunpy.time.parse_time(getHeaderFromFitsFile(filepath)[1]['T_OBS'])

    return (current_time - start_time).seconds

def getDistanceToSunFromSDO(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns distance of the sun to the SDO in meters (DSUN_OBS)
    '''
    
    head = getHeaderFromFitsFile(filepath)[1]
    
    return head["DSUN_OBS"]


def getDataMeanStdDevSkewMedian(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns tuple of (data mean, 
                      standard deviation, 
                      data skew, 
                      data median)
    '''
    
    head = getHeaderFromFitsFile(filepath)[1]
    
    return (head["DATAMEAN"],head["DATARMS"],head["DATASKEW"],head["DATAMEDN"])

def getDataMean(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the data mean from the file header
    '''
    
    return getDataMeanStdDevSkewMedian(filepath)[0]

def getDataRMS(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the data rms from the file header
    '''
    
    return getDataMeanStdDevSkewMedian(filepath)[1]

def getDataSkew(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the data skew from the file header
    '''
    
    return getDataMeanStdDevSkewMedian(filepath)[2]

def getDataMedian(filepath):
    '''
    `filepath` is a path to a FITS file

    Returns the data median from the file header
    '''
    
    return getDataMeanStdDevSkewMedian(filepath)[3]

def getMeanStddevSkewMedian(planet, wavelength):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use
    
    Returns a tuple of (average data mean, 
                        average standard deviation, 
                        average data skew, 
                        average data median)
        for the wavelength
    '''
    total_mean = 0
    total_dev = 0
    total_median = 0
    total_skew = 0
    list_of_files = readListOfFiles(planet, wavelength)
    num_files = len(list_of_files)

    if num_files == 0:
        return (0,0,0,0)
    
    for filepath in list_of_files:
        tup = getDataMeanStdDevSkewMedian(filepath)
        total_mean += tup[0]
        total_dev += tup[1]
        total_skew += tup[2]
        total_median += tup[3]

    return (total_mean/num_files, total_dev/num_files, total_skew/num_files, total_median/num_files)

def getDictOfMeanStddevSkewMedian(planet):
    '''
    `planet` is the planet transit to use

    Returns a dictionary of {
                             wavelength : (average data mean, 
                                           average standard deviation, 
                                           average data skew, 
                                           average data median)
                            }
        for all wavelengths
    '''
    wave_dict = {}
    for wavelength in wavelengths:
        wave_dict[wavelength] = getMeanStddevSkewMedian(planet, wavelength)
        
    return wave_dict

def saveMeanStddevSkewMedianDict(planet):
    '''
    `planet` is the planet transit to use

    Saves a dictionary of {
                           wavelength : (average data mean, 
                                         average standard deviation, 
                                         average data skew, 
                                         average data median)
                          }
        to `personal_data_directory``planet`mean_stddev_skew_median.txt
        for all wavelengths
    '''
    writeToFile(personal_data_directory + planet + "mean_stddev_skew_median.txt", getDictOfMeanStddevSkewMedian(planet))
    
def readMeanStddevSkewMedianDict(planet):
    '''
    `planet` is the planet transit to use

    Returns the saved dictionary of {
                                     wavelength : (average data mean, 
                                                   average standard deviation, 
                                                   average data skew, 
                                                   average data median)
                                    }
        from `personal_data_directory``planet`mean_stddev_skew_median.txt
        for all wavelengths
    '''
    return readFromFile(personal_data_directory + planet + "mean_stddev_skew_median.txt")

mean_stddev_skew_median_dict = {"venus": readMeanStddevSkewMedianDict("venus"), 
                                "mercury": readMeanStddevSkewMedianDict("mercury")}


def readListOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of files from that
    '''
    checkIfInWavelengths(wavelength)
    
    file_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
        if checkTimeBlock(planet, time_block, time):
            file_list.append(filepath)
    return file_list

def readListOfFilesAndTimeSinceStartOfTransit(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of (file, time) from that
    '''
    checkIfInWavelengths(wavelength)
    
    file_and_time_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
        if checkTimeBlock(planet, time_block, time):
            file_and_time_list.append((filepath, time))
    return file_and_time_list

def readExposureTimeOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that
    '''
    
    checkIfInWavelengths(wavelength)
        
    exptime_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
                        
        if checkTimeBlock(planet, time_block, time):
            exptime_list.append(getExpTime(filepath))
        
    return exptime_list

def readActualExposureTimeOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that
    '''
    
    checkIfInWavelengths(wavelength)
        
    exptime_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
                        
        if checkTimeBlock(planet, time_block, time):
            exptime_list.append(getActualExpTime(filepath))
        
    return exptime_list

def readDistancefSunOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the distance to the sun for all the files in the given wavelength
    '''

    checkIfInWavelengths(wavelength)
    distance_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
        
        if checkTimeBlock(planet, time_block, time):
            distance_list.append(getDistanceToSunFromSDO(filepath))

    return distance_list

def readDataMeansOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the data means for all the files in the given wavelength
    '''

    checkIfInWavelengths(wavelength)
    mean_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
        if checkTimeBlock(planet, time_block, time):
            mean_list.append(getDataMean(filepath))
        
    return mean_list

def readTimeSinceStartOfTransitOfFiles(planet, wavelength, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of times from that
    '''
    checkIfInWavelengths(wavelength)
    time_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
        
        if checkTimeBlock(planet, time_block, time):
            time_list.append(time)
           
    return time_list
def getNumberOfPixelsOfFiles(planet, wavelength, split=no_split, rotate=no_rotate, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the number of pixels for all the files in the given wavelength
    '''
    
    checkIfInWavelengths(wavelength)
    
    pixel_count_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
                    
        if checkTimeBlock(planet, time_block, time):
            if split == no_split:
                thingToAppend = pixels
            elif split == inside_split and rotate == transit_rotate:
                thingToAppend = pixels_in_split
            elif split == outside_split and rotate == transit_rotate:
                thingToAppend = pixels_out_of_split
            else:
                thingToAppend = getNumberOfPixels(filepath, split, rotate)
            
            pixel_count_list.append(thingToAppend)
                    
    return pixel_count_list

def readCountsPerSecondOfFiles(planet, wavelength, split=no_split, rotate=no_rotate, time_block=[]):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that
    '''
    
    checkIfInWavelengths(wavelength)
        
    count_list = []
    for filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, 
                   pixels_out_of_split) in readData(planet, wavelength):
                        
        if checkTimeBlock(planet, time_block, time):
            if split == no_split:
                thingToAppend = count
            elif split == inside_split and rotate == transit_rotate:
                thingToAppend = count_in_split
            elif split == outside_split and rotate == transit_rotate:
                thingToAppend = count_out_of_split
            else:
                thingToAppend = getCountsPerSecond(filepath, split, rotate)
            count_list.append(thingToAppend)
        
    return count_list


def findGoodFiles(planet, hour, wavelength):
    '''
    `planet` is the planet transit to use

    `hour` is what hour of the transit to get the data from, from "00" to "23"

    `wavelength` is the wavelength to use
        
    Returns a list of non cutoff files
    '''
    
    if planet == "venus":
        year = "2012"
        month = "06"
        if hour in ["21","22","23"]:
            day = "05"
        else:
            day = "06"
    else:
        year = "2016"
        month = "05"
        day = "09"
    
    checkIfInWavelengths(wavelength)
    
    if wavelength == "cont":
        #Goes to the directory for the given hour and day
        
        os.chdir(hmi_data_directory + year + '/'  + month + '/' + day 
                                                       + '/H' + hour  + '00')
        #Gets a list populated with all the files from the directory
        files = filter(os.path.isfile, os.listdir( os.curdir ) )
        
        #Make sure that there were files on that day
        if files == []:
            print "No files found on day: " + day + " hour: " + hour
            
        good_files = []
        
        for filepath in files:
            good_files.append(hmi_data_directory + year + '/' + month + '/' + day + '/H' + hour + '00/' + filepath)

        return good_files
    else:
        #Goes to the directory for the given hour and day
        os.chdir(aia_data_directory + year + '/' + month + '/' + day + '/H' + hour + '00')
        #Gets a list populated with all the files from the directory
        files = filter(os.path.isfile, os.listdir( os.curdir ) )

        #Looks for all the files that are for the requested wavelength
        correct_wavelength_files = []
        for filepath in files:
            if wavelength + ".fits" in filepath:
                correct_wavelength_files.append(filepath)

        #Make sure that there were files on that day
        if correct_wavelength_files == []:
            print "No files found on day: " + day + " hour: " + hour
            return []

        good_files = []

        #Goes through all the files for the wavelength and puts non-cutoff ones 
        #     in an array to be returned
        for i,filepath in enumerate(correct_wavelength_files):        
            if wavelength == "1600": 
                #All the 1600 wavelength images are cut off, 
                #     so we'll work with all of the cut off ones
                if getHeaderFromFitsFile(filepath)[1]['MISSVALS'] != 0:
                    good_files.append(aia_data_directory + year + '/' + month + '/' + day + '/H' + hour + '00/' + filepath)
            else:
                #Check if the sum of the data array is greater than zero
                #This removes all cutoff images, 
                #     since the black pixels are all equal to -2^15
                if getHeaderFromFitsFile(filepath)[1]['MISSVALS'] == 0:
                    good_files.append(aia_data_directory + year + '/' + month + '/' + day + '/H' + hour + '00/' + filepath)

        return good_files
    
def getAllGoodFiles(planet, wavelength):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    Returns a list of all the non cutoff files
    '''
    checkIfInWavelengths(wavelength)
    
    if planet == "venus":
        good_files =  findGoodFiles("venus","21",wavelength)
        good_files += findGoodFiles("venus","22",wavelength)
        good_files += findGoodFiles("venus","23",wavelength)
        good_files += findGoodFiles("venus","00",wavelength)
        good_files += findGoodFiles("venus","01",wavelength)
        good_files += findGoodFiles("venus","02",wavelength)
        good_files += findGoodFiles("venus","03",wavelength)
        good_files += findGoodFiles("venus","04",wavelength)
        good_files += findGoodFiles("venus","05",wavelength)
    else:
        good_files =  findGoodFiles("mercury","10",wavelength)
        good_files += findGoodFiles("mercury","11",wavelength)
        good_files += findGoodFiles("mercury","12",wavelength)
        good_files += findGoodFiles("mercury","13",wavelength)
        good_files += findGoodFiles("mercury","14",wavelength)
        good_files += findGoodFiles("mercury","15",wavelength)
        good_files += findGoodFiles("mercury","16",wavelength)
        good_files += findGoodFiles("mercury","17",wavelength)
        good_files += findGoodFiles("mercury","18",wavelength)
        good_files += findGoodFiles("mercury","19",wavelength)
        good_files += findGoodFiles("mercury","20",wavelength)
        
    return good_files

def getDictOfAllGoodFilesAllThings(planet, wavelength):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use
    
    Returns a dictionary of             
        filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, pixels_out_of_split)
    '''
    checkIfInWavelengths(wavelength)
    
    good_files = readListOfFiles(planet, wavelength)
    if good_files == []:
        good_files = getAllGoodFiles(planet, wavelength)
            
    d = {}
    
    for filepath in good_files:
        d[filepath] = (getCountsPerSecond(filepath),getTimeSinceStartOfTransit(planet, filepath),
                       getNumberOfPixels (filepath, no_split,      no_rotate),
                       getCountsPerSecond(filepath, inside_split,  transit_rotate),
                       getNumberOfPixels (filepath, inside_split,  transit_rotate),
                       getCountsPerSecond(filepath, outside_split, transit_rotate),
                       getNumberOfPixels (filepath, outside_split, transit_rotate))
    return d

def saveDictData(planet, wavelength, directory=personal_data_directory):
    '''
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `directory` is the directory to save the datafile to
    
    Saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength
        
    The saved data file is named `planet`data`wavelength`.txt

    Returns the dictionary
    '''
    
    checkIfInWavelengths(wavelength)
    
    return writeToFile(directory + planet + 'data' + wavelength + '.txt', 
                       getDictOfAllGoodFilesAllThings(planet, wavelength))
    
def saveAllWavelengthData(planet, directory=personal_data_directory):
    '''
    `planet` is the planet transit to use

    `directory` is the directory to save the datafile to

    For all wavelengths, saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength

        Each saved data file is named `planet`data`wavelength`.txt
    '''
        
    for wavelength in wavelengths:
        saveDictData(planet, wavelength, directory)
        print "Saved data for" + wavelength
            
def saveAllWavelengthDataMultiThreaded(planet):
    '''
    `planet` is the planet transit to use

    For all wavelengths, saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength

        Each saved data file is named `planet`data`wavelength`.txt

    Does this in multiple threads (see section "Multithreading")

    Due to limitations of multithreading, the directory to save data to is the default one from saveDictData
    '''

    #It's so repetitive because putting stuff in a %job does weird stuff with the arguments
    if planet == "venus":
        get_ipython().magic(u'job [saveDictData("venus", "0094")]')
        get_ipython().magic(u'job [saveDictData("venus", "0131")]')
        get_ipython().magic(u'job [saveDictData("venus", "0171")]')
        get_ipython().magic(u'job [saveDictData("venus", "0193")]')
        get_ipython().magic(u'job [saveDictData("venus", "0211")]')
        get_ipython().magic(u'job [saveDictData("venus", "0304")]')
        get_ipython().magic(u'job [saveDictData("venus", "0335")]')
        get_ipython().magic(u'job [saveDictData("venus", "1600")]')
        get_ipython().magic(u'job [saveDictData("venus", "1700")]')
        get_ipython().magic(u'job [saveDictData("venus", "cont")]')
    else:
        get_ipython().magic(u'job [saveDictData("mercury", "0094")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0131")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0171")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0193")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0211")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0304")]')
        get_ipython().magic(u'job [saveDictData("mercury", "0335")]')
        get_ipython().magic(u'job [saveDictData("mercury", "1600")]')
        get_ipython().magic(u'job [saveDictData("mercury", "1700")]')
        get_ipython().magic(u'job [saveDictData("mercury", "cont")]')


def graphLightCurve(planet, wavelength, split=no_split, rotate=no_rotate, popt_one=[], func_one=None, popt_two=[], 
                    func_two=None, 
                    time_block=[], 
                    show_events=False, label="", wavelength_name=True, new_figure=True, scale_to_one=True, shift=0, 
                    scale_to_one_based_on="max", 
                    fontsize=10, numbersize=10, show_graph=True, subplot=None):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)
    
    `func_one` is a function (number, parameters) -> number that is applied to the light curve before graphing it

    `popt_one` is an array of optimized parameters for `func_one`

    `func_two` is a function (number, parameters) -> number that is applied to the light curve before graphing it

    `popt_one` is an array of optimized parameters for `func_two `

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    If `show_events` is True, the graph of the light curve includes lines for various events that happen during the transit
        Shows the ingress and egress with black lines
        Shows midnight (venus transit only) with a green line

    `label` is text to add to the end of the data name that goes in the legend

    If `wavelength_name` is True, the data name that goes in the legend includes the name of the wavelength
    
    If `new_figure` is True, it creates a new matplotlib figure for the light curve

    If `scale_to_one` is True, it scales the light curve to be on a scale of 0 to 1

    `scale_to_one_based_on` decides how to scale the light curve, if `scale_to_one` is True
        "max"         scales so that the maximum element is 1
        "first"       scales so that the first element is 1
        "last"        scales so that the last element is 1
        "not transit" scales so that the average not transit element is 1

    `shift` is added to every datapoint in the array
    
    `fontsize` is the size of the text in the graphed light curve (doesn't include the tick mark labels)
    
    `numbersize` is the size of the tick mark labels in the graphed light curve

    If `show_graph` is True, it graphs the light curve
    
    `subplot` is the argument passed into matplotlib's subplot function. If `subplot == None` then subplot is not used.
    
    Returns a tuple of the array of times and the array of data values that were graphed
    '''
    
    def scaleDataToOutOfOne(planet, data, based_on="max", times=None):
        '''
        `planet` is the planet transit to use

        `data` is the array of data to scale

        Scales the given data to be on a scale of 0 to 1. 

        If `based_on` is not given or is set to "max", it scales to 0 to 1

        If `based_on` is "first", it scales so that the first element is 1

        If `based_on` is "last", it scales so that the last element is 1

        If `based_on` is "not transit", it scales so that the average not transit element is 1
            If this is selected, you must also pass in an array of time elements into `times`
        '''

        #Change to be out of 1
        if based_on=="max":
            max_value = max(data)
        elif based_on=="first":
            max_value = data[0]
        elif based_on=="last":
            max_value = data[-1]
        elif based_on=="not transit":
            if time==None:
                raise Exception("Time parameter needs to be set for scaleDataToOutOfOne with parameter \"not transit\"")
            good_data = getDataFromTimeBlock(planet, ["not transit"], data, times)
            max_value = np.average(good_data)
        else:
            raise Exception("Couldn't scale to out of one because based_on parameter " 
                            + based_on + " is invalid for scaleDataToOutOfOne")
        for i, data_point in enumerate(data):
            data[i] = data_point / max_value

        return data    
    
    checkIfInWavelengths(wavelength)
    
    #Gets all of the sum data points, then divides them by the number of pixels to have averages
    pts = np.divide(readCountsPerSecondOfFiles(planet, wavelength, split, rotate, time_block=time_block),
                    getNumberOfPixelsOfFiles(planet, wavelength, split, rotate, time_block=time_block))
    
    #Gets the list of times
    times = readTimeSinceStartOfTransitOfFiles(planet, wavelength, time_block=time_block)
    
    #Calculate the numbers of sections that make up the split and select either the HMI or AIA split
    if split != None:
        split = split[planet]

        if wavelength == "cont":
            split = split[1]
        else:
            split = split[0]
        split_sections = split[0] * 1.0 / len(split[1]) 
    else:
        split_sections = 1

    if len(popt_one) == 0 and len(popt_two) == 0:
        if scale_to_one:
            pts = scaleDataToOutOfOne(planet, pts, based_on=scale_to_one_based_on, times=times)
            
        for i, data in enumerate(pts):
            #The division needs to be done in the last step.
            #If there is another curve fit, do it there, else do it here.
            pts[i] = pts[i]/split_sections + ((split_sections - 1.0)/ split_sections)
    else:
        #Subtracts the first curve fit from the data
        if len(popt_one) != 0:
            for i, data in enumerate(pts):
                #The division needs to be done in the last step.
                #If there is another curve fit, do it there, else do it here.
                if len(popt_two) != 0:
                    pts[i] = pts[i] - func_one(timeThroughTransit(planet, times[i]),*popt_one) + popt_one[-1]
                else:
                    pts[i] = (pts[i] - func_one(timeThroughTransit(planet, times[i]),*popt_one))/split_sections+popt_one[-1]

        #Subtracts the second curve fit from the data
        if len(popt_two) != 0:
            for i, data in enumerate(pts):
                pts[i] = (pts[i] - func_two(timeThroughTransit(planet, times[i]),*popt_two))/split_sections + popt_two[-1]
        if scale_to_one:
            pts = scaleDataToOutOfOne(planet, pts, based_on=scale_to_one_based_on, times=times)

    #apply shift
    pts = np.add(pts, shift)

    if show_graph:
        if wavelength_name:
            #Name the graph based on the wavelength
            if wavelength == "cont":
                data_name = "HMI Continuum"
            else:
                data_name = "AIA " + wavelength + "Å"  
        else:
            data_name = ""
        data_name += " " + label
        
        #Prepare the canvas
        if new_figure:
            fig = plt.figure()
        else:
            fig = plt.gcf()
        if subplot == None:
            axes = plt.gca()
        else:
            axes = fig.add_subplot(subplot)
            
        #Graph the light curve
        light_curve = LightCurve.create({data_name: pts}, index = times)    
        light_curve.plot(fontsize=numbersize)

        if new_figure:
            plt.xlabel('Time (seconds)', fontsize=fontsize)
            plt.ylabel('Brightness (percent)', fontsize=fontsize)

            if planet == "venus":
                fig.suptitle('Percent Brightness During the 2012 Venus Transit', fontsize=fontsize)
            else:
                fig.suptitle('Percent Brightness During the 2016 Mercury Transit', fontsize=fontsize)

        plt.legend(loc=2,prop={'size':fontsize})

        #Graph the different activities
        if show_events:
            #Ingress
            plt.plot([ingress_start_time[planet], ingress_start_time[planet]], [axes.get_ylim()[0], axes.get_ylim()[1]], 'k-')
            plt.plot([ingress_end_time[planet], ingress_end_time[planet]], [axes.get_ylim()[0], axes.get_ylim()[1]], 'k-')

            #Egress
            plt.plot([egress_start_time[planet], egress_start_time[planet]], [axes.get_ylim()[0], axes.get_ylim()[1]], 'k-')
            plt.plot([egress_end_time[planet], egress_end_time[planet]], [axes.get_ylim()[0], axes.get_ylim()[1]], 'k-')

            if planet == "venus":
                #Midnight
                plt.plot([10800, 10800], [axes.get_ylim()[0], axes.get_ylim()[1]], 'g-')
    
    return (times, pts)

def graphAmountCausedByDistance(planet, wavelength):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    Graphs the light curve and a line representing how much of that was caused by the SDO moving in respect to the Sun
    '''
    
    distance_data = readDistancefSunOfFiles(planet, wavelength)
    times = readTimeSinceStartOfTransitOfFiles(planet, wavelength)
    pts = np.divide(readCountsPerSecondOfFiles(planet, wavelength), getNumberOfPixelsOfFiles(planet, wavelength))
    max_value = max(pts)
    for i, data in enumerate(pts):
        pts[i] = data / max_value
    graphLightCurve(wavelength, scale_to_one=True)
    scaled_data = (np.divide(distance_data, max(distance_data)))
    scaled_data = np.add(scaled_data,max(scaled_data)-min(scaled_data))
    scaled_data = np.multiply(scaled_data, pts[0])
    light_curve = LightCurve.create({"Amount caused by distance": scaled_data}, index = times)    
    light_curve.plot()
    
def graphLightCurveAdjusted(planet, wavelength, show_events=False, use_primary_curve_fit=True, use_secondary_curve_fit=True,
                         label="", wavelength_name=True, new_figure=True, scale_to_one=True, shift=0,
                         force_primary_curve_fit=False, scale_to_one_based_on="max", fontsize=10, numbersize=10, show_graph=True):
    '''
    Parameters that are also in graphLightCurve behave in the same way

    If `use_primary_curve_fit` is True, uses a primary curve fit to adjust the light curve

    If `use_secondary_curve_fit` is True, uses a secondary linear curve fit to adjust the light curve
    
    Returns a tuple of the array of times and the array of data values that were graphed
    '''
    checkIfInWavelengths(wavelength)
    
    if planet == "venus":
        if wavelength in ["1700","1600"] and not force_primary_curve_fit:
            use_primary_curve_fit = False

        if use_primary_curve_fit:
            if wavelength in ["cont"]:
                split_for_primary = outside_split
            else:
                split_for_primary = inside_split

            #Get the times and data for the outside split

            time_block = ["transit"]
            if wavelength == "0304":
                time_block.append("not venus spike")

            times = readTimeSinceStartOfTransitOfFiles(planet, wavelength, time_block=time_block)
            data = readCountsPerSecondOfFiles(planet, wavelength, time_block=time_block,
                                              split=split_for_primary, rotate=transit_rotate)

            num_pixels = getNumberOfPixelsOfFiles(planet, wavelength, time_block=time_block,
                                                  split=split_for_primary, rotate=transit_rotate)
            
            #Adjust the data to be averages instead of totals
            data = np.divide(data, num_pixels)
            
            if wavelength == "0304":
                func_one = quadratic
            else:
                func_one = cubic
                
            #Do curve fitting
            popt_one, pcov = curve_fit(func_one, times, data)

            if split_for_primary == outside_split:
                if wavelength == "cont":
                    num = 1
                else:
                    num = 0
                popt_one = np.divide(popt_one,inside_split[planet][num][0] * 1.0 / len(inside_split[planet][num][1]))
            
        else:
            popt_one=[]
            func_one=None

        if use_secondary_curve_fit:

            #Gets all of the sum data points, then divides them by the number of pixels to have averages
            data = np.divide(readCountsPerSecondOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate),
                             getNumberOfPixelsOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate))

            #Gets the list of times
            times = readTimeSinceStartOfTransitOfFiles(planet, wavelength)

            #Subtracts the curve fit from the data
            if len(popt_one) != 0:
                for i, data_point in enumerate(data):
                    data[i] = data[i] - func_one(timeThroughTransit(planet,times[i]),*popt_one) + popt_one[-1]

            #Get rid of transit points
            new_times = []
            new_data = []

            for i, time in enumerate(times):
                if time < ingress_start_time[planet] or time > egress_end_time[planet]:
                    new_times.append(time)
                    new_data.append(data[i])

            func_two = linear

            #Do curve fitting
            popt_two, pcov = curve_fit(func_two, new_times, new_data)
        else:
            popt_two=[]
            func_two=None

    else:
        if wavelength in ["cont","1700"]:
            if use_primary_curve_fit:
                
                #Gets all of the sum data points, then divides them by the number of pixels to have averages
                data = np.divide(readCountsPerSecondOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate),
                                 getNumberOfPixelsOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate))

                #Gets the list of times
                times = readTimeSinceStartOfTransitOfFiles(planet, wavelength)

                #Get rid of transit points
                new_times = []
                new_data = []

                for i, time in enumerate(times):
                    if checkTimeBlock(planet,["not transit"],time):
                        new_times.append(time)
                        new_data.append(data[i])

                func_one = linear

                #Do curve fitting
                popt_one, pcov = curve_fit(func_one, new_times, new_data)
            else:
                popt_one=[]
                func_one=None
            if use_secondary_curve_fit:
                
                #Gets all of the sum data points, then divides them by the number of pixels to have averages
                data = np.divide(readCountsPerSecondOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate),
                                 getNumberOfPixelsOfFiles(planet, wavelength, split=inside_split, rotate=transit_rotate))

                #Gets the list of times
                times = readTimeSinceStartOfTransitOfFiles(planet, wavelength)

                
                #Subtracts the curve fit from the data
                if len(popt_one) != 0:
                    for i, data_point in enumerate(data):
                        data[i] = data[i] - func_one(timeThroughTransit(planet,times[i]),*popt_one) + popt_one[-1]

                #Get rid of transit points
                new_times = []
                new_data = []

                for i, time in enumerate(times):
                    if checkTimeBlock(planet,["transit"],time):
                        new_times.append(time)
                        new_data.append(data[i])

                func_two = linear

                #Do curve fitting
                popt_two, pcov = curve_fit(func_two, new_times, new_data)
            else:
                popt_two=[]
                func_two=None
        else:
            print "graph adjusted transit doesn't work for mercury " + wavelength
            popt_one=[]
            func_one=None
            popt_two=[]
            func_two=None
    
    #Actually graph the light curve, using the optimal values from the curve fit. 
    #It only uses the proper split
    return graphLightCurve(planet, wavelength, popt_one=popt_one, func_one=func_one, popt_two=popt_two, func_two=func_two, 
                split=inside_split, rotate=transit_rotate, show_events=show_events, label=label,
                wavelength_name=wavelength_name, new_figure=new_figure, scale_to_one=scale_to_one, shift=shift,
                scale_to_one_based_on=scale_to_one_based_on, fontsize=fontsize, numbersize=numbersize, show_graph=show_graph)
    
def graphAllLightCurves(planet, split=no_split, rotate=no_rotate, time_block=[], show_events=False, label="", 
                        wavelength_name=True, 
                        new_figure=True, scale_to_one=True, remove=[], shift_up=False, all_new_figures=False, 
                        scale_to_one_based_on="max", 
                        reverse=False, fontsize=10, numbersize=10, show_graph=True):
    '''
    Parameters that are also in graphLightCurve behave in the same way
    
    If `shift_up` is True, graphs the light curves with 0.0003 gap between them

    If `all_new_figures` is True, graphs each light curve in a new figure, 
        otherwise it graphs all the light curves in a single image

    `remove` is an array of wavelengths not to graph
    '''
    
    for wavelength in remove:
        checkIfInWavelengths(wavelength)

    shift = 0
    
    if reverse:
        r_wavelengths = wavelengths[::-1]
    else:
        r_wavelengths = wavelengths
        
    for wavelength in r_wavelengths:
        if wavelength not in remove:
            if shift_up:
                shift = shift + 0.0003

            graphLightCurve(planet, wavelength, new_figure=new_figure, show_events=show_events,
                            wavelength_name=wavelength_name, label=label,
                            scale_to_one=scale_to_one, time_block=time_block, split=split, rotate=rotate,
                           scale_to_one_based_on=scale_to_one_based_on, fontsize=fontsize, numbersize=numbersize, show_graph=show_graph)
            new_figure=all_new_figures
            show_events=all_new_figures
            

def graphAllLightCurvesAdjusted(planet, use_primary_curve_fit=True, use_secondary_curve_fit=True, show_events=False, 
                                label="", 
                                wavelength_name=True, new_figure=True, scale_to_one=True, remove=[], shift_up=False, 
                                all_new_figures=False,
                                scale_to_one_based_on="max", reverse=False, fontsize=10, numbersize=10, show_graph=True):
    '''
    Parameters that are also in graphLightCurve, graphAllLightCurves, or graphLightCurveAdjusted behave in the same way

    '''
    
    for wavelength in remove:
        checkIfInWavelengths(wavelength)
    
    shift = 0
    
    if reverse:
        r_wavelengths = wavelengths[::-1]
    else:
        r_wavelengths = wavelengths
    
    
    for wavelength in r_wavelengths:
        if wavelength not in remove:
            if shift_up:
                shift = shift + 0.0003
            graphLightCurveAdjusted(planet, wavelength, new_figure=new_figure, show_events=show_events,
                                 use_primary_curve_fit=use_primary_curve_fit, 
                                 use_secondary_curve_fit=use_secondary_curve_fit,
                                 wavelength_name=wavelength_name, label=label,
                                 scale_to_one=scale_to_one, shift=shift,
                                scale_to_one_based_on=scale_to_one_based_on, fontsize=fontsize, numbersize=numbersize, show_graph=show_graph)
            new_figure=all_new_figures
            show_events=all_new_figures
    
def graphWavelengthAndLimbDarkening(planet, wavelength, 
                                    limb_darkening_model="quadratic", limb_darkening_parameters=None, 
                                    depth=None, orbital_period_divider=11.3, semi_major_axis_const=14,
                                    new_figure=True, show_wavelength=True, show_graph=True):
    '''
    Parameters that are also in graphLightCurve behave in the same way

    `depth` increases the magnitude of the predicted light curve change

    `semi_major_axis_const` and `orbital_period_divider` affect the predicted light curve in strange ways

    Graphs the adjusted light curve and a predicted light curve 
        based on the `limb_darkening_model` and `limb_darkening_parameters`

    If `limb_darkening_model` is "quadratic" and `limb_darkening_parameters` is None, 
        it uses parameters from Allen's Astrophysical Quantities
    '''
    
    if planet == "venus":
        if depth==None:
            if wavelength in ["1600","1700"]:
                depth = 1.16
            elif wavelength == "cont":
                depth = 1.45
            else:
                depth = 1.4
        params = batman.TransitParams()
        params.t0 = venus_time_inferior_conjunc      #time of inferior conjunction
        params.per = venus_orbital_period / orbital_period_divider         #orbital period
        params.rp = venus_apparant_radius*depth       #planet radius (in units of stellar radii)
        params.a = semi_major_axis_const #semi-major axis (in units of stellar radii)
        params.inc = 90 - venus_orbital_inclination       #orbital inclination (in degrees)
        params.ecc = venus_eccentricity              #eccentricity
        params.w = venus_longitude_perihelion        #longitude of periastron (in degrees)
        params.limb_dark = limb_darkening_model        #limb darkening model
        if limb_darkening_parameters == None and limb_darkening_model == "quadratic":       #limb darkening coefficients
            #These constants are from Allen's Astrophysical Quantities
            if wavelength == "cont":
                params.u = [0.88, -0.23]
            else:
                params.u = [0.12, 0.33]
        else:
            params.u = limb_darkening_parameters
            
        if params.u == None:
            params.u = []
            
        t = np.linspace(0, 32400, 10000)  #times at which to calculate light curve
        m = batman.TransitModel(params, t)    #initializes model
        flux = m.light_curve(params) #creates light curve
        if new_figure:
            plt.figure()
        plt.plot(t, flux,label='Generated') #plots generated light curve
        if show_wavelength:
            #plots actual light curve
            graphLightCurveAdjusted(planet,wavelength,new_figure=False, scale_to_one_based_on="not transit", show_graph=show_graph) 
    else:
        raise Exception("graphWavelengthAndLimbDarkening not yet implemented for mercury")


def derotateWavelength(planet, wavelength, time_one=0, time_two=30000, file_one=None, file_two=None):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    `time_one` is the time to derotate to

    `time_two` is the time to get an image from to derotate from

    `file_one` is the file to derotate to

    `file_two` is the file to get an image from to derotate from

    Doesn't work for HMI Continuum
    
    Returns the derotated image
    '''
    list_of_maps = []
    if file_one == None:
        list_of_maps.append(sunpy.map.Map(findFileByTimeSinceStartOfTransit(planet, time_one,wavelength)))  
    else:
        list_of_maps.append(sunpy.map.Map(file_one))  
    if file_two == None:
        list_of_maps.append(sunpy.map.Map(findFileByTimeSinceStartOfTransit(planet, time_two,wavelength)))  
    else:
        list_of_maps.append(sunpy.map.Map(file_two))  

    mapcube = sunpy.map.Map(list_of_maps, cube=True)
    derotated = mapcube_solar_derotate(mapcube)
    return derotated.maps[-1]
        
def showDiffBetweenTwoTimes(planet, wavelength, file_one=None, file_two=None, time_one=None, time_two=None, 
                            one_special=False, two_special=False):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    `time_one` is the time to get an image from to difference to

    `time_two` is the time to get an image from to difference from

    `file_one` is the file to get an image from to difference to

    `file_two` is the file to get an image from to difference from

    If `one_special` is True, the array is stored in the first section of the first file

    If `two_special` is True, the array is stored in the first section of the second file
    
    Returns the difference image
    '''
    
    checkIfInWavelengths(wavelength)

    #Gets the two arrays
    if time_one != None:
        file_one = findFileByTimeSinceStartOfTransit(planet, time_one,wavelength)
    if time_two != None:
        file_two = findFileByTimeSinceStartOfTransit(planet, time_two,wavelength)

    if one_special == False:
        image_one = getArrayFromFitsFile(file_one)
    else:
        image_one = np.fliplr(np.flipud(fits.read(file_one)[0][0]))
    if two_special == False:
        image_two = getArrayFromFitsFile(file_two)
    else:
        image_two = np.fliplr(np.flipud(fits.read(file_two)[0][0]))
             
            
    #Finds the differences between the two arrays
    diff_image = np.subtract(image_two, image_one)
    
    #Flips the array so that the orientation is consistent for all wavelengths
    if wavelength == "cont":
        diff_image = np.fliplr(diff_image)
    else:
        diff_image = np.flipud(diff_image)
        
    showArray(diff_image)
    
    return diff_image


def getFourierFrequencyData(planet, wavelength, show_graph=True):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use
    
    If `show_graph` is True, it graphs the frequency data

    Returns the result of a Fourier transform done on the planet and wavelength data
    '''
    
    times, data  = graphLightCurveAdjusted(planet,wavelength,show_graph=False)

    times = np.array(times,np.float64)
    data  = np.array(data, np.float64)
    data  = np.subtract(data,np.average(data))

    seconds_between_points = 24

    # number of signal points
    N = int((max(times) - min(times)) / seconds_between_points)

    x = np.linspace(min(times), max(times), N)
    y = interp1d(times, data)(x)
    yf = fft(y)
    xf = fftfreq(N, seconds_between_points)
    if show_graph:
        xf_n = fftshift(xf)
        yplot = fftshift(yf)
        plt.figure()
        plt.plot(xf_n, 1.0/N * np.abs(yplot))
        plt.grid()
        plt.show()
    return (xf, yf)

def getLowPassedData(planet, wavelength, show_graph=True):
    '''
    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use
    
    If `show_graph` is True, it graphs the low passed data

    Returns the planet and wavelength data with a low pass filter applied
    '''
    
    times, data  = graphLightCurveAdjusted(planet,wavelength,show_graph=False)

    times = np.array(times,np.float64)
    data  = np.array(data, np.float64)
    data  = np.subtract(data,np.average(data))

    seconds_between_points = 24

    # number of signal points
    N = int((max(times) - min(times)) / seconds_between_points)

    x = np.linspace(min(times), max(times), N)
    y = interp1d(times, data)(x)

    fc = 0.03  # Cutoff frequency as a fraction of the sampling rate (in (0, 0.5)).
    b = 0.02  # Transition band, as a fraction of the sampling rate (in (0, 0.5)).
    N = int(np.ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = np.arange(N)

    # Compute sinc filter.
    h = np.sinc(2 * fc * (n - (N - 1) / 2.))

    # Compute Blackman window.
    w = np.blackman(N)

    # Multiply sinc filter with window.
    h = h * w

    # Normalize to get unity gain.
    h = h / np.sum(h)
    y = np.convolve(y, h, mode='same')
    y = np.add(y, data_average)
    if show_graph:
        plt.plot(x,y)
        
    return (x,y)

