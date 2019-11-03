all: OM.so SIM.so
    # Comment here preserves the prior tab
OM.so: OM.cpp
    Rscript --vanilla -e "TMB::compile('OM.cpp','-march=native -O2 -std=gnu++11')"
SIM.so: SIM.cpp
    Rscript --vanilla -e "TMB::compile('SIM.cpp','-march=native -O2 -std=gnu++11')"

clean:
    rm -rf *o 