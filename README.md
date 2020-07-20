**README**

**OVERVIEW**
This is research code implementing "basic" path homology in MATLAB. The key file is pathhomology.m, but this calls reashortestpaths.m. For a more performant and capable implementation building on this one and written by M. Yutin, see a similarly named repository in https://github.com/SteveHuntsmanBAESystems

**CONTENTS**

    pathhomology.m      (path homology)
    reashortestpaths.m  (K shortest paths, using the REA algorithm)
    README-PH.txt       (this file)
    removeloopsisos.m   (useful to preprocess digraphs)

**EXAMPLE** 
To run on an interesting example, try
        
        s = [1,1,2,2,2,3,4,4]; 
        t = [2,4,1,3,4,4,2,3]; 
        D = digraph(s,t); 
        ph = pathhomology(D,4)  % nontrivial homology in dimension 3

**TIPS**
Be careful about the size of input graphs and top dimension. Currently there are no provisions in the code to abort if things will get too big. This is left as an exercise.

Also, note that this code makes use of svds (via rank, though there are commented lines that enable symbolic computation if the toolbox is available), and this sometimes leads to numerical errors, e.g., a negative Betti number. 

**CITE**
If you use this code in your work, please cite it. I would also personally like to hear about your application, evaluation, etc.    

**CONTACT**

    steve.huntsman@baesystems.com

**ACKNOWLEDGEMENTS**
This code is based upon work supported by BAE Systems FAST Labs. 

**LICENSE INFORMATION** 
"Basic" path homology files (including this one) are distributed under the 3-clause BSD license. All licenses are included in the files themselves, and this paragraph is for informational purposes only.

Copyright (c) 2020, BAE Systems. All rights reserved.
