# potts_3d_new
Virtual Cardiac Tissue Model – A Cellular Potts Model for cardiac monolayers that reproduces fibrotic patterns

## Build instructions
`required gcc-9`  
`pip install cython==0.29.32 pyyaml numpy jupyterlab wandb scikit-image feret pandas seaborn numpngw`  
`git clone https://github.com/violinist2802/potts_3d_new.git`  
`cd potts_3d_new`  
`cd VCT`  
`make clean`  
`make`  
`cd ..`  
`python3 setup.py build_ext --inplace`  
open untitled.ipynb  
register on wandb.ai  
copy API key from user_settings in wandb.ai  
login running first cell   
set parameters in ipynb (NCX, NCY....) NCZ = 1  
open sweeps/d20_sweeps/3d.json  
set borders of search for each parameter  
run cell in .ipynb 'with open('sweeps/d20_sweeps/3d.json', 'r')...'  
copy sweep_id from output  
paste in last cell in sweep_id = ''  
run last cell  
look into wandb
