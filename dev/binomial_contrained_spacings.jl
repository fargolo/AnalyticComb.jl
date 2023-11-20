# To do: 
# Substitute l-runs constrained version for conv_bin (DONE)
# Sum up probabilities for inference on SymbolicInference procedure
# Implement weighted case.
# NOTE: Periodicity results in several significant 
### intervals Sharkovskii's theorem

function conv_bin(k,n,d)
    accum = []
    for j in 0:Int(round(n/d))
        if(j>=k-1  || j*d >=n)
            break
        end
        print("\nj is ",j)
        acc = (-1)^j * binomial(k-1,j) * binomial(n-d*j,k)
        print("\naccum is ",acc)
        push!(accum,acc)
    end
    print("\nSum is ")
    sum(accum)
end

counts = map(x->conv_bin(x,10,3),[5,4,3,2])
probs = counts ./ 2^10
test_prob = sum(probs)
plot(probs)