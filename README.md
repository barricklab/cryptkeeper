# cryptkeeper

## Installation

Cryptkeeper is currently in development. Adventurous users can install dependencies on their own through other means. The simplest way to install Cryptkeeper is through conda by:

* Cloning this repository to your machine
* `conda env create -f environment.yml`
* `conda activate cryptkeeper`
* `pip install ./`

It may be necessary to install [rhotermpredict](https://github.com/barricklab/RhoTermPredict), [promotercalculator](https://github.com/barricklab/promoter-calculator), and [ostir](https://github.com/barricklab/ostir) from their repositories to avoid version mismatch as development continues.


## Example

```javascript
cd examples/BPMV
cryptkeeper -i pSMART-LCKan-BPMV1.fna -o pSMART-LCKan-BPMV1
```
