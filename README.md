# APPAC<a href="https://github.com/user-attachments/assets/536458e6-6e1d-4705-a228-8020b54baf3e"/><img src="man/figures/logo.png" align="right" height="96" alte="appac website"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/RuedigerForster/appac/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RuedigerForster/appac/actions/workflows/R-CMD-check.yaml)

[![lint](https://github.com/RuedigerForster/appac/actions/workflows/lint.yaml/badge.svg)](https://github.com/RuedigerForster/appac/actions/workflows/lint.yaml)

<!-- badges: end -->

## Overview

APPAC is a versatile software package designed to improve measurement accuracy. Initially developed for atmospheric pressure peak area correction in gas chromatography, it can also be widely applied to other environmental sensor data. By addressing environmental influences, instrumental bias, and non-linear responses, APPAC provides customized solutions to stabilize and enhance data accuracy across various applications. In addition to outputting correction models, it generates time-based drift factors to effectively mitigate errors caused by baseline shifts, drift, and environmental factors, ensuring accurate and reliable measurements.

## Installation

### Install R

To use this R package from GitHub, the user needs to install **R** at first.  ​**R**​ is the core programming language and environment required to run and install any R package. It must be installed on the user's system.

#### Install R on Windows

1. Go to the [R Project website](https://cran.r-project.org/).
2. Click on ​[Download R for Windows](https://cran.r-project.org/bin/windows/).
3. Click on install R for the first time.
4. Click Download R for Windows. Open the downloaded file.
5. Select the language you would like to use during the installation. Then click OK.
6. Click Next.
7. Select where you would like R to be installed. It will default to your Program Files on your C Drive. Click Next.
8. You can then choose which installation you would like. If your computer is a 64-bit, you can choose the 64-bit User Installation. Then click Next.
9. Then specify if you want to customized your startup or just use the defaults. Then click Next.
10. Then you can choose the folder that you want R to be saved within or the default if the R folder that was created.  Once you have finished, click Next.
11. You can then select additional shortcuts if you would like. Click Next.
12. Click Finish

#### Install R on Linux

1. Go to the [R Project website](https://cran.r-project.org/).
2. Click on ​[Download R for Linux](https://cran.r-project.org/bin/linux/).
3. Click on one of the following links according to your operating system
   * [debian](https://cran.r-project.org/bin/linux/debian/)
   * [fedora](https://cran.r-project.org/bin/linux/fedora/)
   * [redhat](https://cran.r-project.org/bin/linux/redhat/)
   * [suse](https://cran.r-project.org/bin/linux/suse/)
   * [ubuntu](https://cran.r-project.org/bin/linux/ubuntu/)
4. Follow the installation guide to install R.
5. Open R by typing R on your console
6. Run the following commmand to check the installed version of R.
   ```
   version
   ```
   
   

#### Install R on Mac

1. Go to the [R Project website](https://cran.r-project.org/).
2. Click on ​[Download R for macOS]([https://cran.r-project.org/bin/macosx/]
3. Download the latest version of R for macOS and open the `.pkg` file.
4. Follow the installation instructions.
5. Open R by typing R on your console
6. Run the following commmand to check the installed version of R.
   ```
   version
   ```### Install RStudio

While RStudio is not mandatory, it is a popular integrated development environment (IDE) for R. It provides a more user-friendly interface for coding, debugging, and managing R projects. If the user prefers working with a GUI, RStudio can make the process easier.

#### Install RStudio on Windows

1. Visit the [RStudio website](https://posit.co/downloads/).
2. Download the **RStudio Desktop Installer** for Windows.
3. Run the installer and follow the on-screen instructions.

#### Install RStudio on Linux

1. Visit the [RStudio download website](https://posit.co/download/rstudio-desktop/#download).
2. Download the RStudio `.deb` or `.rpm` package according to your operating system and version
3. Install the package

#### Install RStudio on Mac

1. Visit the [RStudio download website](https://posit.co/download/rstudio-desktop/#download).
2. Download the RStudio `.dmg` package according to your MacOS version
3. Install the package


### Install APPAC

The user can install a GitHub package directly in R using tools like `devtools`.
If RStudio is used, the same commands can be run in its console. However, the installation process itself relies only on R.

```r
# Install devtools if not already installed
install.packages("devtools")

# Use devtools to install the package from GitHub
devtools::install_github("RuedigerForster/appac")
```

## Attention :warning:

<p>This package is work in progress and its contents will change frequently. Please
feel free to check it out. What you find here is working and most of it is documented.</p>

<p>Collaborators are welcome.</p>

