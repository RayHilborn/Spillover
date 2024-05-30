README file for Hilborn et al spillover model.

Code has several parts

“Save nine cases.r”  runs the model for different values of the policy, either “effort” or “hr”, and the dispersal parameter and u/umsy parameter and saves each in a “.rdata” file with the name of the policy, dispersal parameter and u/umsy.

It requires “multistep functions v2.r” which contains all the real model code, and “base setup.r” which sets the base parameters.

The code is designed so that multiple steps within a year can be implemented, necessary when dealing with short life history fish, but for this paper there was only one step per year.

After the 18 cases are run then “nine panel plot.r” and “nine panel  results table.r” generate the figures and tables.
