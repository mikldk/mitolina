# mitolina: MITOchondrial LINeage Analysis

An R package (<https://www.r-project.org/>) to perform **MITO**chondrial **LIN**eage **A**nalysis 
by simulating genealogies backwards and 
imposing mtDNA including binary mutations forwards. 
This software operates under the maternal inheritance only model, i.e. that mtDNA is only passed on by mothers to children. 
Numerous analyses are possible, e.g. number of matches and meiotic distance to matches.

See documentation included in package (vignettes and manual) at <https://mikldk.github.io/mitolina/>.

## Installation

You first need `R` (<https://www.r-project.org/>). 
Then you can install `mitolina` from GitHub by using the `remotes` package (<https://CRAN.R-project.org/package=remotes>):

``` r
# install.packages("remotes")
remotes::install_github("mikldk/mitolina", build_opts = c("--no-resave-data", "--no-manual"))
```

The `build_opts` is provided to remove `--no-build-vignettes` such that vignettes are built.

## Getting started

Refer to the included vignettes. You can get an overview of the included vignettes by the following `R` command:

```r
browseVignettes(package = "mitolina")
```

To read a vignette, type:

```r
vignette("introduction", package = "mitolina")
```

## Contribute, issues, and support

Please use the issue tracker at <https://github.com/mikldk/mitolina/issues> 
if you want to notify us of an issue or need support.
If you want to contribute, please either create an issue or make a pull request.


## Disclaimer

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## License

License: GPL-2.

## Build status

Travis CI:

[![Travis-CI Build Status](https://travis-ci.org/mikldk/mitolina.svg?branch=master)](https://travis-ci.org/mikldk/mitolina)

