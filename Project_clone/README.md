This folder contains the project structure that we use to initialize new projects as well as template files for some of the analysis that is common between projects.

If you want to clone this folder to initialize a new project without a need to clone the entire repo do the following:

1. git --version needs to return a number higher than 2.25 if not upgrade git
2. git clone --no-checkout git://github.com/Naturhistoriska/eDNA
3. cd eDNA
4. git sparse-checkout init --cone
5. git sparse-checkout set Project_clone
6. git checkout main
7. cp Project_clone ~/New_projectname

NB! Make sure to update the Project_clone regurarly so you are using the latest templates files.
