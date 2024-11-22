Torch models should be saved with json metadata describing the atom typing scheme and grid resolution/size.  Assuming `model` is a torch model that takes an atom grid as input:

```python
d = {
    'resolution': 0.5,
    'dimension' : 23.5,
    'recmap' : '''AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine Chlorine Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus 
Calcium
Zinc
GenericMetal Boron Manganese Magnesium Iron''',

'ligmap': '''AliphaticCarbonXSHydrophobe 
AliphaticCarbonXSNonHydrophobe 
AromaticCarbonXSHydrophobe 
AromaticCarbonXSNonHydrophobe
Bromine Iodine
Chlorine
Fluorine
Nitrogen NitrogenXSAcceptor 
NitrogenXSDonor NitrogenXSDonorAcceptor
Oxygen OxygenXSAcceptor 
OxygenXSDonorAcceptor OxygenXSDonor
Sulfur SulfurAcceptor
Phosphorus
GenericMetal Boron Manganese Magnesium Zinc Calcium Iron'''
}


extra = {'metadata':json.dumps(d)}
z = torch.zeros((1,28,48,48,48))
        
script = torch.jit.trace(model, z)
script.save(newf,_extra_files=extra)    
```
