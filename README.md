LightCurveCreation by Garrison Taylor

The ipython notebook version offers more functionality 

To get aggregated data, change `personal_data_directory` to the location of the aggregated data directory in this repository. To save new aggregated data, change `aia_data_directory` and `hmi_data_directory` to folders containing year data for those telescopes. These can be downloaded from [JSOC](http://jsoc.stanford.edu/). If using custom `aia_data_directory` and `hmi_data_directory`, you will likely need to save all wavelength data yourself instead of using aggregated data, unless your `aia_data_directory` and `hmi_data_directory` are the same as mine, since the aggregated data includes the absolute paths to those two directories.

== Saving and Reading Aggregated Data ==

- functions -

readData(planet, wavelength)

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    Returns the saved data for that wavelength in list form, sorted by time



writeToFile(filepath, dictionary)

    `filepath` is the filepath to write to

    `dictionary` is the dictionary to save

    Returns the dictionary



readFromFile(filepath)

    `filepath` is the filepath to read from

    Returns the file evaluated

    This will turn, for example, a dictionary in the file 
        into a dictionary in python



== Basic Graphical Tools == 

To change between inline pictures and separate windows for images, change the line starting with %matplotlib to end with either inline, qt, osx, or gtx depending. This only works on the ipython notebook version


- functions -

showImage(filepath)

    `filepath` is a path to a FITS file

    It makes a sunpy map of the given file, and puts the map on screen
    
    Returns the map



showArray(array)

    Graphs the given `array` in greyscale 



showHistogramOfPixelIntensities(array, boxes=50)

    Shows a histogram of elements of the given `array`

    `boxes` is the number of separate bars to sort the elements of the array into
    
    Note: uses lots of ram
        
    Returns the array, with all nan elements set to zero



showSplitOfImageFromFitsFile(filepath, split=no_split, rotate=no_rotate, replace_with_zero=False)

    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    If `replace_with_zero` is true, it will replace the removed elements of the array with zeros, otherwise it graphs a truncated array.

    Graphs the array from the filepath, applying the split and rotate



== Working with FITS Files ==

- functions -

getArrayFromFitsFile(filepath)

    `filepath` is a path to a FITS file

    If the file doesn't exist, wait for it to.

    Returns the array from the given file
    


getHeaderFromFitsFile(filepath)

    `filepath` is a path to a FITS file

    If the file doesn't exist, wait for it to.

    Returns the header from the given file



findFileByTimeSinceStartOfTransit(planet, time, wavelength)

    `planet` is the planet transit to use

    `time` is the time since the start of the transit to get the file from

    `wavelength` is the wavelength to use

    Returns the filepath to the FITS file closest to that time
     

   
getRotatedAndSplitArrayFromFitsFile(filepath, split=no_split, rotate=no_rotate, replace_with_zero=False)

    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    If `replace_with_zero` is true, it will replace the removed elements of the array with zeros, otherwise it returns a truncated array.

    Finds the array from the given file given file, 
        then applies the given split and rotate parameters as needed
        
    Returns the array created



getTotalData(filepath, split=no_split, rotate=no_rotate)

    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)
    
    It ignores the negative items in the array and items with a value 
        greater than the mean + std dev of the wavelength as a whole
    
    Returns the sum of the image data in the fits file



getNumberOfPixels(filepath, split=no_split, rotate=no_rotate)

    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)
    
    Returns the number of non-nan pixels
   

 
getExpTime(filepath)

    `filepath` is a path to a FITS file

    Note: the HMI continuum doesn't have an exposure time, so returns a 1

    Returns the saved exposure time of the file. This may be inaccurate, depending on the processing already done on the file



getActualExpTime(filepath)

    `filepath` is a path to a FITS file

    Returns the actual exposure time of the file, based on the time stamps in the file header



getCountsPerSecond(filepath, split=no_split, rotate=no_rotate)

    `filepath` is a path to a FITS file

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    Returns the number of counts per second from the file, applying the split and rotation angle



getTimeSinceStartOfTransit(filepath)

    `filepath` is a path to a FITS file

    Returns the time in seconds since the start of the recorded transit period



getDistanceToSunFromSDO(filepath)

    `filepath` is a path to a FITS file

    Returns distance of the sun to the SDO in meters (DSUN_OBS)

getPlanetFromFile(filepath)
    
    `filepath` is a path to a FITS file

    Returns either "venus" or "mercury"


getWavelengthFromFile(filepath)
    
    `filepath` is a path to a FITS file



== Data Statistics ==

-- constants --

mean_stddev_skew_median_dict is a dictionary of {"venus": venus_statistics_dictionary, "mercury": mercury_statistics_dictionary}



-- functions --

getDataMeanStdDevSkewMedian(filepath)

    `filepath` is a path to a FITS file

    Returns tuple of (data mean, 
                      standard deviation, 
                      data skew, 
                      data median)



getDataMean(filepath)

    `filepath` is a path to a FITS file

    Returns the data mean from the file header
    


getDataRMS(filepath)

    `filepath` is a path to a FITS file

    Returns the data rms from the file header
    


getDataSkew(filepath)

    `filepath` is a path to a FITS file

    Returns the data skew from the file header
    


getDataMedian(filepath)

    `filepath` is a path to a FITS file

    Returns the data median from the file header
    


getMeanStddevSkewMedian(planet, wavelength)

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use
    
    Returns a tuple of (average data mean, 
                        average standard deviation, 
                        average data skew, 
                        average data median)
        for the wavelength



getDictOfMeanStddevSkewMedian(planet)

    `planet` is the planet transit to use

    Returns a dictionary of {
                             wavelength : (average data mean, 
                                           average standard deviation, 
                                           average data skew, 
                                           average data median)
                            }
        for all wavelengths



saveMeanStddevSkewMedianDict(planet)

    `planet` is the planet transit to use

    Saves a dictionary of {
                           wavelength : (average data mean, 
                                         average standard deviation, 
                                         average data skew, 
                                         average data median)
                          }
        to /home/gtaylor/`planet`mean_stddev_skew_median.txt
        for all wavelengths
    


readMeanStddevSkewMedianDict(planet)

    `planet` is the planet transit to use

    Returns the saved dictionary of {
                                     wavelength : (average data mean, 
                                                   average standard deviation, 
                                                   average data skew, 
                                                   average data median)
                                    }
        from /home/gtaylor/`planet`mean_stddev_skew_median.txt
        for all wavelengths

== Whole Wavelength Lists ==

readListOfFiles(planet, wavelength, time_block=[])
    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of files from that



readListOfFilesAndTimeSinceStartOfTransit(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of (file, time) from that



readExposureTimeOfFiles(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that



readActualExposureTimeOfFiles(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that



readDistancefSunOfFiles(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the distance to the sun for all the files in the given wavelength



readDataMeansOfFiles(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the data means for all the files in the given wavelength



readTimeSinceStartOfTransitOfFiles(planet, wavelength, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of times from that



readCountsPerSecondOfFiles(planet, wavelength, split=no_split, rotate=no_rotate, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Reads the saved data for that wavelength and gets the list of counts per second from that



getNumberOfPixelsOfFiles(planet, wavelength, split=no_split, rotate=no_rotate, time_block=[])

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `split` is the section of the array that is cared about (see section "Split and Rotate")

    `rotate` is the number of degrees to rotate the array by (counter clockwise)

    `time_block` is an array of strings to limit the list to (see section "Time")
    
    Returns a list of the number of pixels for all the files in the given wavelength

== Saving Whole Wavelength Data ==

All functions in this section are very slow, due to the large number of files being processed

-- functions --

findGoodFiles(planet, hour, wavelength)

    `planet` is the planet transit to use

    `hour` is what hour of the transit to get the data from, from "00" to "23"

    `wavelength` is the wavelength to use
        
    Returns a list of non cutoff files



getAllGoodFiles(planet, wavelength)

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    Returns a list of all the non cutoff files



getDictOfAllGoodFilesAllThings(planet, wavelength)

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use
    
    Returns a dictionary of             
        filepath, (count, time, pixels, count_in_split, pixels_in_split, count_out_of_split, pixels_out_of_split)



saveDictData(planet, wavelength, directory='/home/gtaylor/')

    `planet` is the planet transit to use

    `wavelength` is the wavelength to use

    `directory` is the directory to save the datafile to
    
    Saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength
        
    The saved data file is named `planet`data`wavelength`.txt

    Returns the dictionary
    


saveAllWavelengthData(planet, directory='/home/gtaylor/')

    `planet` is the planet transit to use

    `directory` is the directory to save the datafile to

    For all wavelengths, saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength

        Each saved data file is named `planet`data`wavelength`.txt

        
    
saveAllWavelengthDataMultiThreaded(planet)

    `planet` is the planet transit to use

    For all wavelengths, saves a dictionary of 
        all the good files for the wavelength : 
            (their counts per second, their time since start of transit)
        to the data file for that wavelength

        Each saved data file is named `planet`data`wavelength`.txt

    Does this in multiple threads (see section "Multithreading")

    Due to limitations of multithreading, the directory to save data to is the default one from saveDictData

== Graphing Light Curves == 

-- functions --

graphLightCurve(planet, wavelength, split=no_split, rotate=no_rotate, popt_one=[], func_one=None, popt_two=[], func_two=None, time_block=[], show_events=False, label="", wavelength_name=True, new_figure=True, scale_to_one=True, shift=0, scale_to_one_based_on="max")

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

    Graphs the light curve



graphAmountCausedByDistance(planet, wavelength)

    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    Graphs the light curve and a line representing how much of that was caused by the SDO moving in respect to the Sun



graphLightCurveAdjusted(planet, wavelength, show_events=False, use_primary_curve_fit=True, use_secondary_curve_fit=True, label="", wavelength_name=True, new_figure=True, scale_to_one=True, shift=0, force_primary_curve_fit=False, scale_to_one_based_on="max")

    Parameters that are also in graphLightCurve behave in the same way

    If `use_primary_curve_fit` is True, uses a primary curve fit to adjust the light curve

    If `use_secondary_curve_fit` is True, uses a secondary linear curve fit to adjust the light curve
    
    Returns the optimal values for the curve fit
    


graphAllLightCurves(planet, split=no_split, rotate=no_rotate, time_block=[], show_events=False, label="", wavelength_name=True, new_figure=True, scale_to_one=True, remove=[], shift_up=False, all_new_figures=False, scale_to_one_based_on="max")

    Parameters that are also in graphLightCurve behave in the same way
    
    If `shift_up` is True, graphs the light curves with 0.0003 gap between them

    If `all_new_figures` is True, graphs each light curve in a new figure, otherwise it graphs all the light curves in a single image

    `remove` is an array of wavelengths not to graph

graphAllLightCurvesAdjusted(planet, use_primary_curve_fit=True, use_secondary_curve_fit=True, show_events=False, label="", wavelength_name=True, new_figure=True, scale_to_one=True, remove=[], shift_up=False, all_new_figures=False, scale_to_one_based_on="max")
    
    Parameters that are also in graphLightCurve, graphAllLightCurves, or graphLightCurveAdjusted behave in the same way



graphWavelengthAndLimbDarkening(planet, wavelength, limb_darkening_model="quadratic", limb_darkening_parameters=None, depth=None, orbital_period_divider=11.3, semi_major_axis_const=14, new_figure=True)
  
    Parameters that are also in graphLightCurve behave in the same way

    `depth` increases the magnitude of the predicted light curve change

    `semi_major_axis_const` and `orbital_period_divider` affect the predicted light curve in strange ways

    Graphs the adjusted light curve and a predicted light curve based on the `limb_darkening_model` and `limb_darkening_parameters`

    If `limb_darkening_model` is "quadratic" and `limb_darkening_parameters` is None, it uses parameters from Allen's Astrophysical Quantities


linear(x,one,zero)

quadratic(x,two,one,zero)

cubic(x,three,two,one,zero)

quartic(x,four,three,two,one,zero)

quintic(x,five,four,three,two,one,zero)

    Functions that are of the form nth * x^n + n-1th * x^(n-1) + … + one * x + zero

    They work for `func_one` and `func_two` for graphLightCurve



== Difference Images ==

-- functions --

derotateWavelength(planet, wavelength, time_one=0, time_two=30000, file_one=None, file_two=None)

    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    `time_one` is the time to derotate to

    `time_two` is the time to get an image from to derotate from

    `file_one` is the file to derotate to

    `file_two` is the file to get an image from to derotate from

    Doesn't work for HMI Continuum

    Returns the derotated image



showDiffBetweenTwoTimes(planet, wavelength, file_one=None, file_two=None, time_one=None, time_two=None, one_special=False, two_special=False)

    `planet` is the planet transit to use
    
    `wavelength` is the wavelength to use

    `time_one` is the time to get an image from to difference to

    `time_two` is the time to get an image from to difference from

    `file_one` is the file to get an image from to difference to

    `file_two` is the file to get an image from to difference from

    If `one_special` is True, the array is stored in the first section of the first file

    If `two_special` is True, the array is stored in the first section of the second file

    Returns the difference image



== Time ==

`time_block` is an array of strings
    "transit"
    "not transit"
    "ingress"
    "egress"
    "not venus spike"
        don't include stuff during the venus spike

-- constants --

venus_spike_start_time = 9000 

in seconds, the start of spike in 1700Å



venus_spike_end_time = 15000

in seconds, the end of spike in 1700Å



ingress_start_time = {"venus":4020,"mercury":5010} 

in seconds, start of ingress



ingress_end_time = {"venus":5150,"mercury":5280} 

in seconds, end of ingress



egress_start_time = {"venus":26100,"mercury":31440} 

in seconds, start of egress



egress_end_time = {"venus":27350,"mercury":31680} 

in seconds, end of egress



transit_start_time = ingress_end_time



transit_end_time = egress_start_time



-- functions --

timeSinceStartOfTransit(planet, time)

    `planet` is the planet transit to use

    Given a `time` since the start of the recorded period, 
        returns the difference between the given time 
        and when the transit actually started



timeThroughTransit(planet, time)

    `planet` is the planet transit to use

    Given a `time` since the start of the recorded period, 
        returns    the start time if the transit has yet to start,
                   the end time   if the transit has already finished,
                or the current time
        

checkTimeBlock(planet, time_block, time)

    `planet` is the planet transit to use

    `time_block` is an array of strings to limit the list to (see section "Time")
   
    `time` is the time to check

    Checks if `time` is in the `time_block`



getDataFromTimeBlock(planet, time_block, data, times)

    `planet` is the planet transit to use

    `time_block` is an array of strings to limit the list to (see section "Time")

    `data` is an array of generic data points that correspond with the time points in `times`

    `data` and `times` must be of same length    

    Returns the data points that match the `time_block`
    


== Split and Rotate ==

= rotate =

-- constants --

transit_rotate = {"venus":7.4,"mercury":2.0}

    The angle to rotate the picture of the sun so that Venus goes horizontal through the image

no_rotate = 0



= split =

-- constants --

Splits are a dictionary of two tuples of two tuples each that are each in the format (w,[x,y,z,...]). 

That is, it breaks the image into w slices and then uses slices x,y,z,...
The first tuple is for AIA images, the second tuple is for HMI images


inside_split = {"venus":((32,[8]), (16,[3])),"mercury":((256,[159,160,161]), (None))}

    The section of the sun that contains the transit.


outside_split = {"venus":(((32,range(8) + range(9,32)), (16,range(3) + range(4,16)))),"mercury":((256,range(159) + range(162,256)), (None))}

    The section of the sun that doesn't contain the transit


no_split = None



== Wavelengths ==

-- constants --

wavelengths = ["0094","0131","0171","0193","0211","0304","0335","1600","1700","cont"]
    `cont` = HMI intensity continuum, all other are AIA wavelengths



-- functions --

checkIfInWavelengths(wavelength)
    Raises an exception if `wavelength` is not in the list of wavelengths



== Venus Orbit ==

-- constants --

data about Venus

venus_radius = 3760.4 
    
    miles



sun_radius = 432168.6 

    miles



venus_semi_major_axis = 67237909.0

    miles



venus_apparant_radius = 0.023
    
    stellar radii
    ((tangent of angular radius of venus on june 5th 2012)*(1 astronomical unit - radius of sun - radius of earth) / (radius of sun))



venus_longitude_perihelion = 131.53298 
    
    degrees



venus_orbital_period_days = 224.701
    
    days



venus_orbital_period = venus_orbital_period_days * 24 * 60 * 60

    seconds



venus_eccentricity = 0.0067



venus_orbital_inclination = 3.39

    degrees



venus_time_inferior_conjunc = (transit_start_time["venus"] + transit_end_time["venus"] ) / 2

    seconds



== Multithreading ==

This is based on ipython notebook features and thus only works in ipython notebook

-- functions --

%job [stuff to do]

jobs.status() 

    to check on jobs



jobs.traceback(num) 

    for stack trace for dead thread `num`



kill_thread(jobs.all[num]) 

    to kill thread `num`



for thread in jobs.running:
    kill_thread(thread)


    Kill all running threads


