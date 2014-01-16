This package contains the best 4 solutions submitted into the "Collective Minds Exploration Challenge" competition.

You can find the competition task problem statement here: http://community.topcoder.com/longcontest/?module=ViewProblemStatement&rd=15765&pm=12621.
It is necessary to read the problem statement in order to understand what the solutions are doing.

The solutions are located in files place1.cpp, place2.cpp, place3.cpp and place4.cpp.

In order to use each of them, it needs to be compiled first (examples below always use place1 solution):

g++ --std=c++0x -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3 -Wl,--stack,4194304 -o place1 place1.cpp

The obtained executable can be run on the example test case using the official competition's visualizer.
You can read more about it here: http://community.topcoder.com/contest/problem/StructureRecognition/manual.html.
For the goals of this package it was slightly modified.

Visualizer is to be run as follows:

java -jar tester.jar -folder data -training example_train.csv -test example_test.csv -exec place1 -vis result

The solution's total score will be displayed as the last line of the visualizer's output.
The visualizations will be saved to the "result" folder.
