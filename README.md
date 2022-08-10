# vital_rscripts_git
r scripts used for the data analysis of the vital project




# Notes to start working with github on linux 

1 create repository online
2 copy repository to linux computer (using a personal access token for authorization) (git clone https://user_name:<token>@github.com/user_name/repository.git) see: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token

git clone https://d-schlaeppi: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@github.com/d-schlaeppi/vital_rscripts_git.git

3 access folder in terminal (cd .. to go back one folder)
cd Documents/vital_rscripts_git

4 At the beginning of each session update from the online repository
git pull

5 When finished to safe it all online
git add *   (says that all the files shall be collected for updating/uploading)
git commit (commit to the changes)
in the terminal text editor explain your changes and then exit (control X), yes for safe, and Enter to confirm name
git push https://github.com/d-schlaeppi/vital_rscripts_git.git

