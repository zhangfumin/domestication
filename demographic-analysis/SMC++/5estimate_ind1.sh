#date_200827_ind1_estimate
#wangmeixia

export PATH=/build/Cellar/anaconda2/bin:$PATH
source activate /build/Cellar/anaconda2/
conda activate smcpp

#take the ind1 population for example
smc++ estimate --cores 2 --timepoints 10  20000 -o model_time 1e-8 ind_chr*_input
smc++ estimate --cores 2 --timepoints 10  20000 -o model_knots --knots 15 1e-8 ind_chr*_input
smc++ estimate --cores 2 -o model 1e-8 ind_chr*_input                 

