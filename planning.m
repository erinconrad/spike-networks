%{

Ideas:
- what's up with HUP105-106: it looks like mscohere is very sensitive to
any dc offset of the signals
       - I can just remove the offset, BUT I SHOULD CONFIRM MSCOHERE IS
       CORRECT
- recollect HUP075 1-35?
- need at least 20 patients
- try delta band? https://www.biorxiv.org/content/biorxiv/early/2020/06/01/2020.05.27.118802.full.pdf
 this paper found pre-IED increased oscillatory power in the delta and
 theta bands

- decide if I should remove bad channels
    - I would either need to do this at the outset, applying it to all
    spikes, or I would need to figure out how to combine adjacency matrices
    across spikes when I have different subsets of channels for each spike

Questions:


%}