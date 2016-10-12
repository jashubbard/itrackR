itrackR Basics
================
Jason Hubbard



itrackR is an R package for high-level analyses of eyetracking data. For now it is only compatible with EDF files from SR-Research Eyelink eyetrackers. **This is currently a work-in-progress.** Some functions are missing documentation, but you can see an example of a full analysis with `vignette('itrackR')`. You will first need to install SR Research's API, found here: (https://www.sr-support.com/forumdisplay.php?17-EyeLink-Display-Software) and the [edfR](http://github.com/jashubbard/edfR) package in order to import the edf files. These aren't on CRAN yet, so everything should be installed via [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```r
install.packages('devtools')
devtools::install_github('jashubbard/itrackR')
```
--- 


###NEWS (10-12-2016)
- Version 0.9
- Some documentation now included. 
- Totally revamped code for handling pupil data (using a sqlite backend)
- Total rewrite of epoching, aggregating, and plotting fixations

---- 


This will acquaint you with some of the basics of itrackR. There are a number of function for doing high-level analyses of eyetracking data.

Example Data
------------

Some sample data is included with itrackR. This includes 2 EDF files from an SR-Research Eyetracker, and a behavioral file (saved as a .rds). You can load them using `itrackR.data`:

``` r
library(itrackR)
datapath <- itrackr.data('path')
edf_files <- itrackr.data('edfs')
beh <- itrackr.data('beh')
```

`datapath` will point to the data folder wherever itrackR is installed:

``` r
datapath
```

    ## [1] "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/itrackR/extdata/"

`edfs` is a list of the 2 edfs found in that folder:

``` r
edf_files
```

    ## [1] "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/itrackR/extdata//104_exp.edf"
    ## [2] "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/itrackR/extdata//105_exp.edf"

`beh` is a data frame of the behavioral data for the same 2 subjects.

|   ID|  Block|  Trial|  Task|  Conflict|  Targetpos|  Distractorpos|
|----:|------:|------:|-----:|---------:|----------:|--------------:|
|  104|      1|      1|     1|         1|          2|             11|
|  104|      1|      2|     1|         0|          4|             NA|
|  104|      1|      3|     1|         0|          2|             NA|
|  104|      1|      4|     1|         0|          6|             NA|
|  104|      1|      5|     1|         0|          4|             NA|
|  104|      1|      6|     1|         0|          2|             NA|
|  104|      1|      7|     1|         1|          3|              8|
|  104|      1|      8|     1|         1|          6|              9|
|  104|      1|      9|     1|         1|          6|              9|
|  104|      1|     10|     1|         1|          3|              7|

Loading Data
------------

We start by initializing the itrackR object and loading the data. This can be done in a couple different ways. If we have a list of all edf files (in the current working directory, or with full path names) we can do it like this.

``` r
z <- itrackr(edfs=edf_files)
```

    ## Loading file /Library/Frameworks/R.framework/Versions/3.2/Resources/library/itrackR/extdata//104_exp.edf
    ## Loading Events (59201 events across 1325 trials)....
    ## Done with 104_exp
    ## 
    ## Loading file /Library/Frameworks/R.framework/Versions/3.2/Resources/library/itrackR/extdata//105_exp.edf
    ## Loading Events (61362 events across 1325 trials)....
    ## Done with 105_exp
    ## 
    ## 
    ## Done processing all files!

``` r
#Alternatively, we can provide the path and a search pattern to find all edfs in a certain folder:
#z <- itrackr(path=datapath, pattern='*.edf')
```

Object Structure
----------------

The `itrackR` object consists of fields for each relevant event. Each one is a data frame. The `ID` fields specifies each subject. The subject ID is formed from extracting only the numeric data from the EDF file.

`z$fixations`

|     |  eyetrial|  sttime|  entime|   gavx|   gavy|   ID|  fixation\_key|
|-----|---------:|-------:|-------:|------:|------:|----:|--------------:|
| 6   |         1|       8|     159|  493.9|  433.5|  104|              1|
| 10  |         1|     172|    2817|  528.0|  409.5|  104|              2|
| 16  |         1|    3036|    3106|  520.4|  396.8|  104|              3|
| 20  |         1|    3112|    4423|  527.8|  399.7|  104|              4|
| 26  |         1|    4713|    8625|  525.2|  401.3|  104|              5|
| 32  |         1|    8956|   10101|  534.7|  406.3|  104|              6|
| 38  |         1|   11171|   11327|  500.1|  437.3|  104|              7|
| 42  |         1|   11349|   11813|  529.6|  406.7|  104|              8|
| 48  |         1|   12732|   12788|  527.7|  416.9|  104|              9|
| 52  |         1|   12799|   13577|  532.1|  411.5|  104|             10|

`z$saccades`

|     |  eyetrial|  sttime|  entime|   gstx|   gsty|   genx|   geny|   avel|    pvel|   ID|  saccade\_key|
|-----|---------:|-------:|-------:|------:|------:|------:|------:|------:|-------:|----:|-------------:|
| 8   |         1|     160|     171|  497.2|  430.9|  512.7|  417.5|   50.7|    64.9|  104|             1|
| 14  |         1|    2818|    3035|  532.9|  403.6|  522.5|  360.7|  295.0|  1024.5|  104|             2|
| 18  |         1|    3107|    3111|  522.3|  404.2|  527.7|  400.7|   36.7|    39.0|  104|             3|
| 24  |         1|    4424|    4712|  522.2|  391.9|  509.8|  375.7|  273.1|  1001.0|  104|             4|
| 30  |         1|    8626|    8955|  523.4|  400.6|  520.7|  384.3|  236.9|   928.0|  104|             5|
| 36  |         1|   10102|   11170|  536.3|  399.6|  574.7|  542.5|  242.9|   775.6|  104|             6|
| 40  |         1|   11328|   11348|  496.9|  444.5|  514.9|  408.7|   57.6|    96.3|  104|             7|
| 46  |         1|   11814|   12731|  534.1|  401.9|  507.9|  730.0|  264.4|   830.8|  104|             8|
| 50  |         1|   12789|   12798|  528.1|  415.5|  521.0|  399.4|   52.0|    61.4|  104|             9|
| 56  |         1|   13578|   14196|  535.1|  409.2|  543.4|  391.1|  202.2|   965.2|  104|            10|

`z$messages` shows the messages that were sent to Eyelink during the experiment

|     |  eyetrial|  sttime| message                    |   ID|  message\_key|
|-----|---------:|-------:|:---------------------------|----:|-------------:|
| 3   |         1|       1| !MODE RECORD CR 1000 2 1 L |  104|             1|
| 5   |         1|      55| BASELINE\_START BLOCK 1    |  104|             2|
| 166 |         1|   60053| BASELINE\_END BLOCK 1      |  104|             3|
| 172 |         2|   85998| !MODE RECORD CR 1000 2 1 L |  104|             4|
| 174 |         2|   86036| TRIALSTART                 |  104|             5|
| 175 |         2|   86036| BLOCK 1                    |  104|             6|
| 176 |         2|   86036| TRIAL 1                    |  104|             7|
| 181 |         2|   87098| STIMONSET                  |  104|             8|
| 198 |         2|   89710| RESPONSE                   |  104|             9|
| 204 |         3|   89755| !MODE RECORD CR 1000 2 1 L |  104|            10|

ROIs
----

Much of the analysis depends on specifying regions of interest (ROIs). We can then determine if fixations lie within these ROIs. First we specify all possible ROIs that may occur in an experiment. The function `radialCoords` makes it easy to specify a set of evenly-spaced coordinates arranged in a ring. We will create elliptical ROIs in this example. `roiFlower` makes it easy to rotate the ellipses to make a flower-like pattern.

``` r
#generate coordinates for our ROIs
innercoords <- radialCoords(x=512,y=384,numpoints=6,radius=240);
outercoords <- radialCoords(512, 384,6, 280,starting_angle=30); #larger radius, starting w/ 30 degree offset

#specify rotations of ellipses
angles <- roiFlower(12)
```

We use `makeROIs` to specify them. We can make elliptical or circular ROIs. First we make the inner ones. We can check our progress using `plot_rois`. Note the plots have the origin at the upper-left. This means that ROI \#1 is actually at the 12 o'clock position if we viewed the plot in the regular orientation.

``` r
#make elliptical ROIs and plot them
z <- makeROIs(z,innercoords,shape='ellipse',xradius=60, yradius=120, angles=angles[c(1,3,5,7,9,11)])
plot_rois(z)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-10-1.png)

Now we add the outer ROIs. We make sure and specify the `append` option, and also provide names. For now, names are limited to numbers. If you don't specify them, it will just use 1...n.

``` r
#make elliptical ROIs and plot them
z <- makeROIs(z,outercoords,shape='ellipse',xradius=60, yradius=120, angles=angles[c(2,4,6,8,10,12)], names=7:12, append=T)
plot_rois(z)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-11-1.png)

Finally let's include a central, circular ROI.

``` r
#coordinates have to be a matrix:
centercoords <- matrix(c(512,384),nrow=1)

z <- makeROIs(z,centercoords,shapes='circle',radius=65, names=13, append=T)
plot_rois(z)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-12-1.png)

Plotting
--------

Once the ROIs are added, we can easily make scatterplots of fixations for each subject. This allows us to find calibration issues.

``` r
plot(z, zoom=TRUE)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-13-1.png)

You can also plot everyone's data on one plot using the `oneplot` argument.

``` r
plot(z,zoom=TRUE,oneplot=TRUE)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-14-1.png)

Merging with behavioral data
----------------------------

Ideally, we send a message to Eyelink on every trial in order to identify it in the EDF file. Every time the eyetracker is started and stopped, we call this a separate trial. In this example, every trial the message "BLOCK X" and "TRIAL X" was sent to eyelink. We can see the trial-wise information in the header of our object:

`z$header`

``` r
knitr::kable(head(z$header, 10))
```

|  eyetrial|  starttime|  endtime|  duration|   ID|  first\_sample|
|---------:|----------:|--------:|---------:|----:|--------------:|
|         1|          1|    60056|     60055|  104|         220490|
|         2|      85998|    89719|      3721|  104|         220490|
|         3|      89755|    92002|      2247|  104|         220490|
|         4|      92041|    93849|      1808|  104|         220490|
|         5|      93886|    95800|      1914|  104|         220490|
|         6|      95834|    98145|      2311|  104|         220490|
|         7|      98185|   100203|      2018|  104|         220490|
|         8|     100234|   103344|      3110|  104|         220490|
|         9|     103374|   105712|      2338|  104|         220490|
|        10|     105738|   107800|      2062|  104|         220490|

Next we specify index variables that uniquely identify trials. This should be present in both the edf and behavioral file. `set_index` searches through the messages, finds the relevant ones, and extracts the numeric data. `set_index` can take a regular expression to find anything that matches this pattern, and `numeric.only` tells it to ignore any text (e.g., "BLOCK "). The variable names are stored in `z$indexvars` and the information is added to `z$header`.

`find_messages` is for pulling other message information from the EDF file. Here we want to have the timestamps of the stimulus onset and the response so we can refer to them later.

`add_behdata` merges the behavioral file with the eye data, based on the index variables. This only works if we run `set_index` first (so we have `Block` and `Trial` variables in the behavioral data frame, and in our itrackR object).

``` r
#find messages to use as our index variables (to merge with our behavioral data)
z<- set_index(z,varnames=c('Block','Trial'), patterns=c('^BLOCK [0-9]*','^TRIAL [0-9]*'), numeric.only=T)

#find messages that specify the onset of events, extract the timestamps
z <- find_messages(z,varnames=c('STIMONSET','RESPONSE'), patterns=c('STIMONSET','RESPONSE'), timestamp=T)
#merge with behavioral data
z <- add_behdata(z,beh,append=F)
```

Now `z$beh` contains the behavioral data that matches with the eyetracking data, based on the index variables you specified (`Block` and `Trial`). It adds the variable `eyetrial`, which is also found in the `header`,`fixations`, `saccades`, `blinks`, and `messages` data frames in the itrackR object. It also adds the timestamps that you requested using `find_messages`:

|   ID|  Block|  Trial|  Task|  Conflict|  Targetpos|  Distractorpos|  eyetrial|  starttime|
|----:|------:|------:|-----:|---------:|----------:|--------------:|---------:|----------:|
|  104|      1|      1|     1|         1|          2|             11|         2|      85998|
|  104|      1|     10|     1|         1|          3|              7|        11|     107836|
|  104|      1|     11|     1|         0|          2|             NA|        12|     110025|
|  104|      1|     12|     1|         1|          4|             11|        13|     111913|
|  104|      1|     13|     1|         1|          6|             12|        14|     114633|
|  104|      1|     14|     1|         1|          4|              9|        15|     117898|
|  104|      1|     15|     1|         1|          6|              8|        16|     120763|
|  104|      1|     16|     1|         1|          5|             11|        17|     123646|
|  104|      1|     17|     1|         0|          2|             NA|        18|     126331|
|  104|      1|     18|     1|         0|          4|             NA|        19|     128505|

We can also see that the header has been updated

`z$header`

``` r
knitr::kable(head(z$header, 10))
```

|  eyetrial|  starttime|  endtime|  duration|   ID|  first\_sample|  Block|  Trial|  STIMONSET|  RESPONSE|
|---------:|----------:|--------:|---------:|----:|--------------:|------:|------:|----------:|---------:|
|         1|          1|    60056|     60055|  104|         220490|     NA|     NA|         NA|        NA|
|         2|      85998|    89719|      3721|  104|         220490|      1|      1|      87098|     89710|
|         3|      89755|    92002|      2247|  104|         220490|      1|      2|      90825|     91998|
|         4|      92041|    93849|      1808|  104|         220490|      1|      3|      93110|     93846|
|         5|      93886|    95800|      1914|  104|         220490|      1|      4|      94966|     95798|
|         6|      95834|    98145|      2311|  104|         220490|      1|      5|      96903|     98142|
|         7|      98185|   100203|      2018|  104|         220490|      1|      6|      99255|    100197|
|         8|     100234|   103344|      3110|  104|         220490|      1|      7|     101299|    103341|
|         9|     103374|   105712|      2338|  104|         220490|      1|      8|     104438|    105709|
|        10|     105738|   107800|      2062|  104|         220490|      1|      9|     106802|    107797|

Now we can plot a subset of our data based on the behaivoral data. Just use the `condition` argument when plotting. Note it uses a similar syntax as `subset`. Here we're plotting data only from the first 5 blocks:

``` r
plot(z,zoom=TRUE,oneplot=TRUE, condition = Block<=5)
```

    ## Warning: Removed 303 rows containing missing values (geom_point).

![](itrackR_files/figure-markdown_github/unnamed-chunk-19-1.png)

Drift Correction
----------------

It looks like subject 104 is off-center, due to poor calibration. We can correct for this using `drift_correct`. We can optionally specify a grouping variable (from the behavioral data) so that correction is done separately for each level of that variable. In this case, let's perform correction for each subject and block. The threshold specifies the minimum amount of movement detected before we actually do any correction.

``` r
z <- drift_correct(z,vars='Block',threshold=15)
plot(z,zoom=T)
```

![](itrackR_files/figure-markdown_github/unnamed-chunk-20-1.png)

Much better!

Determining Fixation/Saccade "Hits"
-----------------------------------

Next we code whether each fixation and saccade "hit" any of the ROIs using `calcHits`

``` r
z <- calcHits(z)
```

Note that the `fixations` data frame now has binary vectors for each ROI specifying whether the fixation hit that item or not:

``` r
knitr::kable(head(z$fixations, 10))
```

|  fixation\_key|   ID|  eyetrial|  sttime|  entime|   gavx|   gavy|  roi\_1|  roi\_2|  roi\_3|  roi\_4|  roi\_5|  roi\_6|  roi\_7|  roi\_8|  roi\_9|  roi\_10|  roi\_11|  roi\_12|  roi\_13|
|--------------:|----:|---------:|-------:|-------:|------:|------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|--------:|--------:|--------:|--------:|
|              1|  104|         1|       8|     159|  493.9|  410.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              2|  104|         1|     172|    2817|  528.0|  386.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              3|  104|         1|    3036|    3106|  520.4|  374.0|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              4|  104|         1|    3112|    4423|  527.8|  376.9|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              5|  104|         1|    4713|    8625|  525.2|  378.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              6|  104|         1|    8956|   10101|  534.7|  383.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              7|  104|         1|   11171|   11327|  500.1|  414.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              8|  104|         1|   11349|   11813|  529.6|  383.9|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|              9|  104|         1|   12732|   12788|  527.7|  394.1|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|
|             10|  104|         1|   12799|   13577|  532.1|  388.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|

This is not terribly useful if your task-relevant ROI changes positions on each trial. You can use `mapROIs` map your experiment-wide ROIs (1,2,3...13) to trial-specific ROIs ('target','distractor'). You just need a variable in your behavioral data that specifies the number of the relevant ROI. Here, `Targetpos` specifies the target location, and `Distractorpos` specifies the distractor location:

``` r
z <- mapROIs(z,names=c('target','distractor'),indicators=c('Targetpos','Distractorpos'))
```

Now we can see a `target_hit` and `distractor_hit` variable in our fixation data frame

``` r
knitr::kable(head(z$fixations, 10))
```

|  fixation\_key|   ID|  eyetrial|  sttime|  entime|   gavx|   gavy|  roi\_1|  roi\_2|  roi\_3|  roi\_4|  roi\_5|  roi\_6|  roi\_7|  roi\_8|  roi\_9|  roi\_10|  roi\_11|  roi\_12|  roi\_13|  Targetpos|  Distractorpos|  target\_hit|  distractor\_hit|
|--------------:|----:|---------:|-------:|-------:|------:|------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|--------:|--------:|--------:|--------:|----------:|--------------:|------------:|----------------:|
|              1|  104|         1|       8|     159|  493.9|  410.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              2|  104|         1|     172|    2817|  528.0|  386.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              3|  104|         1|    3036|    3106|  520.4|  374.0|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              4|  104|         1|    3112|    4423|  527.8|  376.9|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              5|  104|         1|    4713|    8625|  525.2|  378.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              6|  104|         1|    8956|   10101|  534.7|  383.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              7|  104|         1|   11171|   11327|  500.1|  414.5|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              8|  104|         1|   11349|   11813|  529.6|  383.9|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|              9|  104|         1|   12732|   12788|  527.7|  394.1|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|
|             10|  104|         1|   12799|   13577|  532.1|  388.7|       0|       0|       0|       0|       0|       0|       0|       0|       0|        0|        0|        0|        1|         NA|             NA|           NA|               NA|

Saving data
-----------

Next we probably want to do statistics on our eyetracking data. We want to have our behavioral data merged with the eye data, including our ROI "hits". Just use `eyemerge` to pull out the relevant information. Any eyetracking data that does not match behavioral data will still be included (all behavioral variables will just be `NA`).

``` r
fixes <- eyemerge(z,'fixations')

#including only some behavioral variables. ID and indexvars are always included
saccs <- eyemerge(z,'saccades',behdata=c('Task'))

#by default only mapped ROIs are included. Here we can include all 13 rois, plus the mapped ones
fixes_all <- eyemerge(z,'fixations',all.rois = T)
```

We can also subset based on our behavioral variables using the `condition` argument:

``` r
fixes_first5 <- eyemerge(z,'fixations',condition = Block <= 5)
max(fixes_first5$Block,na.rm=T)
```

    ## [1] 5

Timeseries plots
----------------

Sometimes you want to see the tendency of the eyes to look at a particular ROI over time, relative to some event. To look at this, we first determine epochs around our event of interest, in this case, `STIMONSET`. We first run `epoch_fixations` for each ROI that we're interested in.

``` r
#start at stimulus onset, going 700ms after that point. Bin the data into 25ms time bins.
#repeate for the target and the distractor ROIs we created
z <- epoch_fixations(z,c('target','distractor'),event='STIMONSET',start = 0, end = 700, binwidth = 25)
```

Next we want to visualize these timeseries using `plot_fixation_epochs`. We can generate separate lines for each level of a factor (specified in the behavioral data) or for different ROIs. You can specify variables that define the different lines, as well as the rows and columns in separate panels. Behind the scenes, the function first aggregates based on `ID` and the factors your specify, then across subjects.

``` r
#plot the timeseries data for fixations to target, separately for the Conflict and Task conditions
plot_fixation_epochs(z,event='STIMONSET',rois=c('target'),group=c('Conflict'),cols='Task')
```

    ## $plot

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](itrackR_files/figure-markdown_github/unnamed-chunk-28-1.png)

    ## 
    ## $data
    ## Source: local data frame [116 x 7]
    ## Groups: Conflict, Task, timepoint [116]
    ## 
    ##    Conflict  Task timepoint    roi   val grouping  color
    ##       <int> <int>     <dbl>  <chr> <dbl>   <fctr> <fctr>
    ## 1         0     1         0 target     0        0      0
    ## 2         0     1        25 target     0        0      0
    ## 3         0     1        50 target     0        0      0
    ## 4         0     1        75 target     0        0      0
    ## 5         0     1       100 target     0        0      0
    ## 6         0     1       125 target     0        0      0
    ## 7         0     1       150 target     0        0      0
    ## 8         0     1       175 target     0        0      0
    ## 9         0     1       200 target     0        0      0
    ## 10        0     1       225 target     0        0      0
    ## # ... with 106 more rows

``` r
#plot fixations to target and disctractor for the same conditions
#you must specify 'roi' as one of the plotting variables for it to work.
plot_fixation_epochs(z,event='STIMONSET',rois=c('target','distractor'),group=c('roi','Conflict'),color='Conflict',cols='Task')
```

    ## $plot

    ## Warning: Removed 60 rows containing missing values (geom_point).

![](itrackR_files/figure-markdown_github/unnamed-chunk-28-2.png)

    ## 
    ## $data
    ## Source: local data frame [232 x 7]
    ## Groups: roi, Conflict, Task [8]
    ## 
    ##           roi Conflict  Task timepoint   val     grouping  color
    ##         <chr>    <int> <int>     <dbl> <dbl>       <fctr> <fctr>
    ## 1  distractor        0     1         0   NaN distractor.0      0
    ## 2  distractor        0     1        25   NaN distractor.0      0
    ## 3  distractor        0     1        50   NaN distractor.0      0
    ## 4  distractor        0     1        75   NaN distractor.0      0
    ## 5  distractor        0     1       100   NaN distractor.0      0
    ## 6  distractor        0     1       125   NaN distractor.0      0
    ## 7  distractor        0     1       150   NaN distractor.0      0
    ## 8  distractor        0     1       175   NaN distractor.0      0
    ## 9  distractor        0     1       200   NaN distractor.0      0
    ## 10 distractor        0     1       225   NaN distractor.0      0
    ## # ... with 222 more rows

``` r
#plot difference waves (target - distractor fixations). Plot on separate rows insetead of columns.
#This example doesn't make much sense because distractors aren't present on no-conflict trials
plot_fixation_epochs(z,event='STIMONSET',rois=c('target','distractor'),group=c('roi','Conflict'),rows='Task',type='difference')
```

    ## $plot

    ## Warning: Removed 59 rows containing missing values (geom_point).

![](itrackR_files/figure-markdown_github/unnamed-chunk-28-3.png)

    ## 
    ## $data
    ## Source: local data frame [116 x 7]
    ## Groups: roi, Conflict, Task [4]
    ## 
    ##           roi Conflict  Task timepoint   val     grouping        color
    ##         <chr>    <int> <int>     <dbl> <dbl>       <fctr>       <fctr>
    ## 1  difference        0     1         0   NaN difference.0 difference.0
    ## 2  difference        0     1        25   NaN difference.0 difference.0
    ## 3  difference        0     1        50   NaN difference.0 difference.0
    ## 4  difference        0     1        75   NaN difference.0 difference.0
    ## 5  difference        0     1       100   NaN difference.0 difference.0
    ## 6  difference        0     1       125   NaN difference.0 difference.0
    ## 7  difference        0     1       150   NaN difference.0 difference.0
    ## 8  difference        0     1       175   NaN difference.0 difference.0
    ## 9  difference        0     1       200   NaN difference.0 difference.0
    ## 10 difference        0     1       225   NaN difference.0 difference.0
    ## # ... with 106 more rows

``` r
#Repeat, but on a subset of the data. Say only when Conflict==1
plot_fixation_epochs(z,event='STIMONSET',rois=c('target','distractor'),group=c('roi','Conflict'),rows='Task',type='difference',condition = Conflict==1)
```

    ## $plot

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](itrackR_files/figure-markdown_github/unnamed-chunk-28-4.png)

    ## 
    ## $data
    ## Source: local data frame [58 x 7]
    ## Groups: roi, Conflict, Task [2]
    ## 
    ##           roi Conflict  Task timepoint         val     grouping
    ##         <chr>    <int> <int>     <dbl>       <dbl>       <fctr>
    ## 1  difference        1     1         0  0.25000000 difference.1
    ## 2  difference        1     1        25  0.05555556 difference.1
    ## 3  difference        1     1        50  0.03846154 difference.1
    ## 4  difference        1     1        75  0.02777778 difference.1
    ## 5  difference        1     1       100  0.02500000 difference.1
    ## 6  difference        1     1       125  0.02380952 difference.1
    ## 7  difference        1     1       150  0.02380952 difference.1
    ## 8  difference        1     1       175 -0.04545455 difference.1
    ## 9  difference        1     1       200 -0.36708861 difference.1
    ## 10 difference        1     1       225 -0.46642857 difference.1
    ## # ... with 48 more rows, and 1 more variables: color <fctr>
