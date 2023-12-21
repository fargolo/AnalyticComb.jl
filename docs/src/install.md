
# Installation  

Python package sympy is required for some functionalities. SymPy.jl is reexported in `AnalyticComb.jl`.   

```
$python -m pip install --upgrade pip
$pip install sympy
```

Then, from Julia:  
```
julia> # type the right bracket to enter pkg REPL ']'
pkg>add AnalyticComb
```

Or

```
julia>using Pkg; Pkg.add("AnalyticComb")
```

