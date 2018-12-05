# Kaggle Santa 2018 - TRP

## MIP
to run MIP. Start Julia
```
using TRP
```

you probably first need to add it by: `] dev TRP` .

Then you can run:

```
 TRP.main_mip("combined","mip", 0,200; N=100)
```

This reads the submission file "combined.csv" and if there is an improvement it will be stored in `submissions/mip.csv`.
Starting MIP from 0 up to 200 with a path length of 100. Current stride is 45 but will be customizable.