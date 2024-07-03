# HTP Beta regression example

1. To create a new Project by cloning this repository from Github: 

      - Open RStudio and go to File > New Project... > Version Control > Git.  
      - Enter the repository URL.  
      - You may want to modify Project directory name and/or location.  
      - Click on Create Project.  
      - You may need to enter Github Username and PAT (or run `usethis::create_github_token()`).  
      - Once project opens in RStudio, in R console run `renv::init()` to initialize project and install required packages  
      (renv should already be installed: `install.packages("renv")`).  

2. To open an existing Project:  

      a. Navigate to and double-click .Rproj file.  
         or     
      b. Open RStudio and go to File > Open Project... > Select .Rproj file.  

3. Once Project is open in RStudio, open R script(s) and work through analysis steps in script. 
