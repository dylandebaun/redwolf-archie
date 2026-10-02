ms term meaning placeholderms
ms {{N_TOTAL}} total haplotypes = N_focal + 1 + N_ref derived
-I 3 {{N_FOCAL}} 1 {{N_REF}}       HYBRID / RWEXTINCT(always 1) / COYOTE      from sample lists
-en 0 1...     present HYBRID size     SIZE_HYBRID
-en {{T_RWEXTINCT}} 2...        RW collapse to extinction     T_RWEXTINCT,SIZE_RWEXTINCT
-en 0 3 1             COYOTE = reference size (fixed at 1)          none
-es {{T_ADMIX}} 1...       admixture: frac stays in HYBRID, rest → ghost     T_ADMIX, ADMIX_STAY
-ej {{T_ADMIX}} 4 3        ghost joins COYOTE (= introgression source)        shares T_ADMIX
-ej {{T_RWMERGE}} 2 1        RWEXTINCT merges into HYBRID          T_RWMERGE
-en {{T_ANC}} 1..           ancestral size change                 T_ANC, SIZE_ANC
-ej {{T_ANC}} 1 3            HYBRID merges into COYOTE (root)shares      T_ANC
