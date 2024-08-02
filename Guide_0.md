# Guide 0:

## A. Download required documents

1. Download the workshop's datasets by clicking [here](https://github.com/TriLab-bioinf/WORKSHOP2024_3/archive/refs/heads/main.zip).
2. Uncompress the downloaded `WORKSHOP2024_3-main.zip` file. You will see a new directory named `WORKSHOP2024_3-main`.
3. Open Rstudio
4. Go to `File > New Project > Existing Directory > Browse` and select the `WORKSHOP2024_3-main` folder and click on `Create Project` button.
5. On the console run the following R commands to set the environment with all required R packages:

```
renv::activate()

renv::restore(confirm = FALSE)

renv::snapshot()
```
