//readme.txt
//HW 4 GUI
//Author: Yuding Ai
//Penn ID: 31295008
//My OS: Ubuntu 16.04
//Date: 2017.03.28

Implemented Features:
1. All the features from Code Requirements: 
    -Main Menu
    -File Dialog for Open and Save files
    -Toolbars for open and save
    -Create 3 Dockable panes:
        -Camera Dock for manully adjust camera settings, using Spinboxes
        -Rasterize Dock for rasterize obj with options, using QRadioButtons
        -Filter Dock for filtering img with options, using QCheckboxes
        (to save space, Rasterize dock and filter dock are tablified together
         so that it will only display one of the two in the front at a time)
    -Camera Dock takes care of modifying the camera
    -Create 2 Buttons for the Rasterize Dock and Filter Dock respectively,
     and when the button is clicked, it will rasterize the opened obj or filtered
     the opened ppm img

2. Finishing the first two extra credit requirements:
    -Allow the user to open an image (rather than obj)and filter it and ignoring
     the camera and rasterization .
    -Implemented the interactive camera using KeypressEvent such that:
        - Rotate left: left arrow Key
        - Rotate right: right arrow Key
        - Increase FOV: 1
        - Decrease FOV: 2
        - Translate right: D
        - Translate left: A
        - Translate up: Q
        - Translate down: E
        - Zoom in: W
        - Zoom out: S

Thanks very much!
