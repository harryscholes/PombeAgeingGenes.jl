"""
    AbstractFile

A type representing a specific type of file to be loaded.

Subtypes must be constructed using `@file F`

```
@file FooFile
load(F::File{FooFile}) = ...
const Foo = FooFile("path", "file")
```
"""
abstract type AbstractFile end

abstract type AbstractFileCollection <: AbstractFile end

filepath(F::AbstractFile) = F.filepath

Base.isfile(F::AbstractFile) = isfile(filepath(F))

Base.open(F::AbstractFile) = open(filepath(F))

"""
    @file F

Creates a composite type `T` with the name `<F>File` and a constant `F`. Useful constructors and methods are also implemented.
"""
macro file(F, args...)
    T = Symbol("$(F)File")

    esc(quote
        struct $T <: AbstractFile
            filepath::String
        end

        const $F = ($T)(eval.($args)...)
    end)
end
