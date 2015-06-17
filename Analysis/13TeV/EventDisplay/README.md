Event Display code
=== 

This is a ROOT macro to make event display that shows both reco and gen information
for ttbar and signal events. Before running the macro, make two directories where 
the generated event display files will be saved. 
```
mkdir ttbar 
mkdir sig
```

The macro has some drawing options which can be turned on/off as an argument to the 
following function,
```
void EventDisplayForBaby(bool sig=false, bool truth=true, bool fatjet=true, char* Region="SR0")
```

This is what these arguements do :   

```
sig     : draw signal or ttbar (sig=false draws ttbar)
truth   : draw truth info or not (truth=true draws truth information) 
fatjet  : draw fatjet or not (fatjet=true draws truth information) 
Region  : Selection defined in ../PassSelection.h  
```

To run the code, do 
```
root -b -q EventDisplayForBaby.C++\(false,true,true,\"SR0\"\)
```

It will generate the plots and print out information about fat jets and skinny jets.  

#### Possible improvements 


