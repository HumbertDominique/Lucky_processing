# Lucky_processing
Software use to track, process and stack star pictures. The main goal is to process focused and defocussed images that were taken on the same sensor.

## Context
This piece of software was developped in the scope of a bachelor thesis in microengineering at HEIG-VD (heig-vd.ch) formaly named "Caractérisation du télescope de l'observatoire solaire IRSOL (Tessin)" (Characterization of the IRISOL solar observatory telescope (Ticino))

## About tracking 
- speckle tracking is implemented
- centroïd tracking is not implemented
## Limitations
- .fits only
- Processing is done all at once. No batches to fit into the RAM.
## Status
Working, could be improved

## Quickstart

The procedure to run Lucky Processing is as follows:
-  either:
(a) Place the files and folders as shown on section 5.3.2
-  install Python
-  install required packages in the work folder in Python with
pip install -r requierments.txt
-  change the filename and file path in flat_proc.py and functions.py
-  either:
   -  Preffered:
   -  input ```flat_proc.py```
   -  input python ```main_focus.py```
   -  input python ```main_defocus.py```
-  or:
   -  run python in the work directory
   -  input ```import flat_proc``` (this will process the flat frame)
   -  input ```import functions as fnc```
   -  input ```fnc.lucky_process_focus()``` (this will run the process on focused images)
   -  input ```fnc.lucky_process_defocus()``` (this will run the process on defocused images)
   -  
It is recommended to only run the program on a few frames and display the output in order to dial in the parameters. It is possible to start at any point of the process as long as the previous
steps were already performed.

Display can be coded inside ```Lucky_process``` or using ```compare.py```

## Examples
### From the command line

```
flat_proc.py
    python main_focus.py 1 10       # Correct the ROI in functions.py
    python main_focus.py r = 0.01   # process the whole batch and selects the 1% best frames
    python main_defocus.py          # process the whole batch
```

### From the command line in python
```
flat_proc.py
    py # Start python
    import flat_proc
    import functions as fnc
    fnc.lucky_process_focus(1, 10)
    exit()                              # Correct the ROI in functions.py
    py # Start python
    import functions as fnc
    fnc.lucky_process_focus(1, 10)
    fnc.lucky_process_focus(r = 0.01)   # process the whole batch and selects the 1% best frames
    fnc.lucky_process_defocus()         # process the whole batch
    ```