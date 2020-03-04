# LDsim
An R/shiny simulation app for teaching linkage disequilibrium.

For less than 10 users a web version is available at 
https://cjbattey.shinyapps.io/LDsim/ .  

For larger classes I recommend running locally by pasting the following code into an R console. 

```
 pkgs <- c("plyr","ggplot2","shiny")
 dl_pkgs <- subset(pkgs,!pkgs %in% rownames(installed.packages()))
 if(length(dl_pkgs)!=0){
   for(i in dl_pkgs) install.packages(i)
 }
 library(shiny)
 runGitHub(username="cjbattey",repo="LDsim")
```

