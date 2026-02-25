# Description of the WOMBATmid ocean biogeochemical model
## Subroutine - "update_from_bottom"

```
        (\___/)  .-.   .-. .--. .-..-..---.  .--. .-----.
        / o o \  : :.-.: :: ,. :: '' :: .; :: .; :'-. .-'
       (   "   ) : :: :: :: :: :: .. ::   .':    :  : :
        \__ __/  : '' '' ;: :; :: :; :: .; :: :: :  : :
                  '.,'.,' '.__.':_;:_;:___.':_;:_;  :_;

 World Ocean Model of Biogeochemistry And Trophic-dynamics (WOMBAT)
```

_Contact Pearse J. Buchanan for any questions_

_Pearse.Buchanan@csiro.au_

---

The subroutine `generic_WOMBATmid_update_from_bottom` moves sinking organic material from the water column into the sediment pools.
It is at this point that the model performs permanent burial of sinking organic matter if desired.

---

### Permanent burial of particulates.

If `do_burial = .true.`, we compute the fraction of incident sinking particualte carbon, iron, silicon and $CaCO_3$ that is permanently buried in the sediments. This permanently buried fraction is effectively removed from the model and therefore is not accumulated within the sedimentary pools.

The fraction buried is calculated according to Equation 3 of [Dunne et al. (2007)](https://doi.org/10.1029/2006GB002907):

$F_{bury} = 0.013 \cdot 0.53 \dfrac{(f_{org})^{2}}{\left(7 + f_{org}\right)^{2}}$

where $f_{org}$ is the rain rate of organic carbon detritus on the seafloor in [mmol C m<sup>-2</sup> s<sup>-1</sup>].

As organic matter rains down at a more rapid rate, the fraction of incident organic carbon, organic iron, organic silicon and $CaCO_3$ that is buried increases.
