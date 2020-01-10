"""
    @in(df, column, collection)

Return a `DataFrame` of the rows in dataframe `df` whose `column` value are in
`collection`.
```
df = DataFrame(A=1:10, B="a")
@in(df, :A, [1,10])
```
"""
macro in(df, column, collection)
    esc(quote
        $df[map(x->in(x, $collection), $df[:,$column]),:]
    end)
end

function counter(xs)
    d = Dict{eltype(xs),Int}()

    for x = xs
        haskey(d, x) ? d[x] += 1 : d[x] = 1
    end

    return d
end
