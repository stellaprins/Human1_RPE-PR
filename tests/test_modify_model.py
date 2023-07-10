import pytest
import os
from pathlib import Path
import sys
sys.path.insert(0,'.')
os.chdir('.') 
from src.modify_model import add_all_EX_rxns, remove_compartment, fix_compartment_dict, add_id_suffix, close_PR_EX, close_EX, create_permutation_dicts, change_bounds
from cobra.io import read_sbml_model

@pytest.fixture
def basis_model(scope='module'):
    return read_sbml_model(Path().cwd() / "tests/data/Human-GEM.xml")

@pytest.fixture
def model(scope='module'):
    return read_sbml_model(Path().cwd() / "tests/data/model_06_MeanExpression_RPE__VoigtEtAl2019_.xml")

# add_all_EX_rxns()
def test_length_add_all_EX_rxns(basis_model, model):
    new_model = add_all_EX_rxns(basis_model, model)
    EX_rxns_basis_model = [r for r in basis_model.reactions if len(r.products)==0]
    EX_rxns_new_model = [r for r in new_model.reactions if len(r.products)==0]
    assert len(EX_rxns_basis_model) == len(EX_rxns_new_model)
def test_id_set_add_all_EX_rxns(basis_model, model):
    new_model = add_all_EX_rxns(basis_model, model)
    EX_rxns_basis_model = [r.id for r in basis_model.reactions if len(r.products)==0]
    EX_rxns_new_model = [r.id for r in new_model.reactions if len(r.products)==0]
    assert set(EX_rxns_basis_model) == set(EX_rxns_new_model)

# remove_compartment()    
def test_remove_compartment(model):
    compartment = list(model.compartments.keys())[0]
    new_model = remove_compartment(model.copy(), compartment)
    assert  new_model.compartments != model.compartments
def test_compartment_number_post_remove_compartment(model):
    compartment = list(model.compartments.keys())[0]
    new_model = remove_compartment(model.copy(), compartment)
    assert len(new_model.compartments.keys()) == len(model.compartments.keys())-1
    
# fix_compartment_dict()    
def test_fix_compartment_dict(model):
    # add values to dictionary keys (necessary for valid sbml model)
    new_model = fix_compartment_dict(model)
    tf = '' in new_model.compartments.values()
    assert tf == False
    
# add_id_suffix()
def test_reactions_add_id_suffix(model):
    add_id_suffix(model,'_suffix')
    assert set([r.id.endswith('_suffix') for r in model.reactions]) == {True}
def test_metabolites_add_id_suffix(model):
    add_id_suffix(model,'_suffix')
    assert set([m.id.endswith('_suffix') for m in model.metabolites]) == {True}
    
# close_PR_EX()  
def test_close_PR_EX_no_PR(model):
    # open all reactions
    for r in model.reactions:
        r.bounds=(-1000,1000)
    # close PR_EX reactions
    model = close_PR_EX(model)
    EX_PR_bounds_zero = set([r.bounds for r in model.reactions if '_PR' in r.id if len(r.products) == 0]) == {(0,0)}
    EX_PR_bounds_nonexistent = set([r.bounds for r in model.reactions if '_PR' in r.id if len(r.products) == 0]) == set()
    assert EX_PR_bounds_zero or EX_PR_bounds_nonexistent == True
def test_close_PR_EX_PR(model):
    # open all reactions
    for r in model.reactions:
        r.bounds=(-1000,1000)
    # close PR_EX reactions
    add_id_suffix(model,'_PR')
    model = close_PR_EX(model)
    EX_PR_bounds_zero = set([r.bounds for r in model.reactions if '_PR' in r.id if len(r.products) == 0]) == {(0,0)}
    EX_PR_bounds_nonexistent = set([r.bounds for r in model.reactions if '_PR' in r.id if len(r.products) == 0]) == set()
    assert EX_PR_bounds_zero or EX_PR_bounds_nonexistent == True
    
# close_EX()  
def test_close_EX(model):
    # open all reactions
    for r in model.reactions:
        r.bounds=(-1000,1000)
    # close EX reactions
    model = close_EX(model)
    assert set([r.bounds for r in model.reactions if len(r.products) == 0]) == {(0,0)}

def test_create_permutation_dicts_N():
    test_dict = {'met1':[1,2], 'met2':[1,2,3,4], 'met3':[1,2,3,4,5,6]}
    test_perm_dict = create_permutation_dicts(test_dict)
    assert len(test_perm_dict) == 2*4*6
    
def test_change_bounds_new_bounds(model):
    model_new = change_bounds(model,{model.reactions[0].id:(-666,666)})
    assert model_new.reactions[0].bounds == (-666,666)