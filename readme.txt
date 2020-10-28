# Quaternion Based Image Shadow Removal

This repository contains the code for the paper "[Quaternion Based Image Shadow Removal]"


Software - Matlab R2019b
----------------------------------------
Usage

code/stroke.m - to get user input (strokes on lit and shadow regions)

code/main.m   - for shadow removal

code/color_transfer2 - Code to transfer color from one image to another [1].
code/qrot3d.c        - Code for 3d quaternion rotation [2].

code/imgs/original  - input shadow image
code/imgs/strokes   - input user strokes on shadow and non-shadow regions
code/imgs/detection - shadow detection input, if available 
code/imgs/removal   - shadow removal output


** all images are in png format
----------------------------------------

References:

[1] Reinhard, E.; Adhikhmin, M.; Gooch, B.; Shirley, P. Color transfer between images. IEEE Computer graphics and applications, 21(5):34-41, 2001.

[2] Steven Michael (2018). qrot3d -- Quaternion Rotation (https://www.mathworks.com/matlabcentral/fileexchange/7107-qrot3d-quaternion-rotation), MATLAB Central File Exchange.