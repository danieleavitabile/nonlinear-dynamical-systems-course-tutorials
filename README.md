[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7330169.svg)](https://doi.org/10.5281/zenodo.7330169)

# Nonlinear Dynamical Systems Part 3 - Tutorial

This repository contains tutorials for the course Nonlinear Dynamical Systems, taught at Vrije Universiteit Amsterdam.

The tutorials cover Part 3 of the course, which gives elements of pattern formation in parabolic PDEs. The course is aimed to 3rd year mathematics students who have followed an introductory courses on Dynamical Systems and Numerical Analysis, have little background on PDEs, and no background in Functional Analysis. Hence its emphasis is on parallels/contrasts with ODEs, formal analysis, and numerics.

The tutorials, which are aimed for a 3 4-week course, build slowly and cover the basics of numerical bifurcation analysis and pattern formation. Students are expected to self-study, which is why solutions and code are distributed upfront. 

## For the impatients
* To have a peek at the tutorial's content, open `Tutorial1/tutorial-01.pdf` and `Tutorial1/solutions-01-code.html`, and similar.
* A more advanced tutorial, aimed at postgraduate students, can be found [here](https://zenodo.org/record/3821169#.Y3X94S8w0gw) 

## Getting started
1. Go to the `Tutorial1` folder, which contains links to files in the `src` folder (no need to enter the latter folder yet). 
1. Open `Tutorial1/tutorial-01.pdf` and start working: make a brand new folder on your computer for Tutorial 1, and populate it with your own code (keep solutions on a side).
1. Solutions to Tutorial 1 are available to you in an html file, so you can click on `Tutorial1/solutions-01-code.html` and browse at once code and generated graphs. Some questions are answered with pen and paper, and solutions are in `Tutorial1/solutions-01-analysis.pdf`.
1. You should be able to progress even if you get stuck on a question: try to get some tips from the html file.
1. Some tasks are longer than others, and involve substantial amount of coding, which is another reasons for providing solutions. Even if you find it too cumbersome to write code for every single question, you are expected to at least understand well the solutions before progressing.
1. The code that generated `Tutorial1/solutions-01-code.html` is in `src/Tutorial1/Solutions/Code`, and you may need access to it for some tasks. While you are welcome to reuse functions in this directory in other tutorials, my strong advice is to write your own files first, and then amend your functions after you have seen the solutions, if you find it useful. Remember that there is not a unique way to accomplish a task.

You are welcome to send me pull requests if you correct errors in the tutorial questions or solutions: all sources, including the LaTeX files for the questions are in the `src/` directory, so you can modify them too (you'll be acknowledged on the tutorials).

I am keen to include Julia and Python versions of the solutions, so do not hesitate to contact me if you want to help.

## Citation
If you use ideas or code from these tutorials, please cite them as follows

```
@software{Tutorials_for_Nonl_Dyn_Sys_2022,
  author = {Avitabile, Daniele},
  doi = {10.5281/zenodo.7330169},
  month = {11},
  title = {{Tutorials for Nonlinear Dynamical Systems - Part 3}},
  url = {https://github.com/danieleavitabile/nonlinear-dynamical-systems-course-tutorials},
  version = {0.0.1},
  year = {2022}
}
```
