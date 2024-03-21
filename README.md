# ecRPE-PR models
Adding enzyme-constraints to genome-scale metabolic models (GEMs) of the outer retina 

## Background

### Goal
Creation of a GEMs of the outer retina (i.e. the RPE and PR layers) that can be used to apply constraint-based reconstruction and analysis (COBRA) methods to investigate hypotheses about eye metabolism. 

## Creation of cell-specific GEMs
We created cell-specific models using an implementation of Cost Optimization Reaction Dependency Assessment (CORDA2; Schultz et al., 2017; Schultz & Qutub, 2016) for Python implemented by Christian Diener (https://github.com/resendislab/corda).

**Template model**

We used Human1 (Robinson et al., 2020) as template model (28-06-2023; https://github.com/SysBioChalmers/Human-GEM).

**Single-cell expression data:**
* RPE: Voigt et al., 2019 (Supplementary data, File pnas.1914143116.sd01.xlsx, column **mean_expr_7**)
* PR rod: Lukowski et al., 2019 (Supplementary data, Dataset EV1, File EMBJ-38-e100811-s003.xlsx, **Average Rod PR C0, C1, C2, C3, C4, C7**)
* PR cone: Lukowski et al., 2019 (Supplementary data, Dataset EV1, File EMBJ-38-e100811-s003.xlsx, **C10 Cone PR**)

1. **Discretize expression into the following expression confidence levels:**
   - &nbsp;&nbsp;**-1** &nbsp;&nbsp;for not present 			    (< .00001)
   - &nbsp;&nbsp;**0** &nbsp;&nbsp;for unknown confidence 		(missing value)
   - &nbsp;&nbsp;**1** &nbsp;&nbsp;for low confidence 		    (.00001 > Q25)
   -  &nbsp;&nbsp;**2**&nbsp;&nbsp; for medium confidence 		 (Q25 > Q75)
   -  &nbsp;&nbsp;**3** &nbsp;&nbsp;for high confidence 		   (> Q75)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: for Lukowski PR rod and cone data most values fell into expression confidence levels **-1**, **2**, or **3**.

2. **Mapping of expression confidence levels to reaction confidence levels**
   -  Using the *reaction_confidence* function from the CORDA package; "and" is evaluated by the minimum confidence and "or" by the maximum confidence.

3. **Reconstruction of cell-specific models**
   - objective function: biomass equation (so CORDA removes reactions in pathways that are not contributing to the production of biomass)
   - exchange bounds unchanged from downloaded Human1 model

## Creation of the RPE-PR Models in COBRApy

1. **Load PR and RPE Models.**
2. **Add Back Lost Exchange Reactions and ATP Hydrolysis Reaction.**
3. **Give Metabolite, Reaction, and Compartment IDs Unique Suffixes**
   - *_RPE* or *_PR* depending on the cell-type.
4. **Create RPE-PR Interface for the PR Model**
   - Add RPE-PR interface reactions by duplicating reactions involving the extracellular space (*[e_PR]*) and assigning new suffixes:
     - Metabolite ID suffix: *e_RPE_PR*
     - Reaction ID suffix: *_PR-RPE*
     - Compartment: *[e_PR] -> [e_RPE_PR]*
5. **Create RPE-PR Interface for the RPE Model**
   - Add RPE-PR interface reactions by duplicating reactions involving the extracellular space (*[e_RPE]*) and assigning new suffixes:
     - Metabolite ID suffix: *e_RPE_PR*
     - Reaction ID suffix: *_RPE-PR*
     - Compartment: *[e_RPE] -> [e_RPE_PR]*
6. **Fuse RPE and PR Models**
   ```python
   # fuse models
   mod_RPE_PR = mod_RPE.copy()
   mod_PR_copy = mod_PR.copy()
   mod_RPE_PR.id = mod_RPE.id + '_' + mod_PR.id
   mod_RPE_PR.name = mod_RPE.id + '_' + mod_PR.id
   mod_RPE_PR.add_reactions(mod_PR_copy.reactions)

7. **Delete duplicated reactions**
   - Delete one duplicate of duplicated interface reactions (that solely take place in the [e_RPE_PR] compartment), and change reaction IDs of remaining reaction to *_eRPE_PR* 

8. **Fix compartment dictionary**
   - Make sure all compartment dictionary keys are associated to values (otherwise SBML model is not valid). 

9. **Make sure added reactions are added to their 'subsystem' groups**
```ruby
    # add all reactions to their subsystem group
    g_dict = {g.name:g.id for g in model.groups}
    for r in model.reactions:
        model.groups.get_by_any(g_dict[r.subsystem])[0].add_members([r])
```

10. **Save model in SBML format**

## A schematic overview

![image](https://github.com/stellaprins/ecRPE-PR/assets/30465823/b8a92360-1acf-4b3b-8a0f-a75ed0c04a36)

## References
Diener, C. (2023, September 1). CORDA for Python. https://github.com/resendislab/corda

Lukowski, S. W., Lo, C. Y., Sharov, A. A., Nguyen, Q., Fang, L., Hung, S. S., Zhu, L., Zhang, T., Grünert, U., Nguyen, T., Senabouth, A., Jabbari, J. S., Welby, E., Sowden, J. C., Waugh, H. S., Mackey, A., Pollock, G., Lamb, T. D., Wang, P., … Wong, R. C. (2019). A single-cell transcriptome atlas of the adult human retina. The EMBO Journal, 38(18). https://doi.org/10.15252/EMBJ.2018100811

Robinson, J. L., Kocabaş, P., Wang, H., Cholley, P. E., Cook, D., Nilsson, A., Anton, M., Ferreira, R., Domenzain, I., Billa, V., Limeta, A., Hedin, A., Gustafsson, J., Kerkhoven, E. J., Svensson, L. T., Palsson, B. O., Mardinoglu, A., Hansson, L., Uhlén, M., & Nielsen, J. (2020). An atlas of human metabolism. Science Signaling, 13(624). https://doi.org/10.1126/SCISIGNAL.AAZ1482

Schultz, A., Mehta, S., Hu, C. W., Hoff, F. W., Horton, T. M., Kornblau, S. M., & Qutub, A. A. (2017). Identifying cancer specific metabolic signatures using constraint-based models. Pacific Symposium on Biocomputing, 0, 485–496. https://doi.org/10.1142/9789813207813_0045

Schultz, A., & Qutub, A. A. (2016). Reconstruction of Tissue-Specific Metabolic Networks Using CORDA. PLOS Computational Biology, 12(3), e1004808. https://doi.org/10.1371/JOURNAL.PCBI.1004808

Voigt, A. P., Mulfaul, K., Mullin, N. K., Flamme-Wiese, M. J., Giacalone, J. C., Stone, E. M., Tucker, B. A., Scheetz, T. E., & Mullins, R. F. (2019). Single-cell transcriptomics of the human retinal pigment epithelium and choroid in health and macular degeneration. Proceedings of the National Academy of Sciences of the United States of America, 116(48), 24100–24107. https://doi.org/10.1073/PNAS.1914143116

## Contact
Phil: p.luthert@ucl.ac.uk (Institute of Ophthalmology, University College London, London, UK)

Stella: s.prins@ucl.ac.uk (Institute of Ophthalmology, University College London, London, UK)


