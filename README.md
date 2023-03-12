# Real-Time Raytracing Renderer

This is a group project done as part of the COnputer Graphics course at TU Delft. It is a RT renderer that has support for reflections, hard/soft shadows, texture mapping and many others.

## My Contribution

I was mainly responsible for the BVH optimization. First I implemented a more basic version which was followed by the SAH binning version. Using the BVH, the render times improved by a significant margin (speed up based on the amount of triangles in the scene)

## Build and Run

To run this project, open it as a folder in Visual Studio, wait for everything to setup and then click the start button at the top.
You can play around with the setting in the GUI. By default, you are in Rasterized mode so you will have to change that if you want to see the real-time raytracing