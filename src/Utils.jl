struct ConvergenceError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::ConvergenceError)
    print(io, "ConvergenceError: ")
    print(io, e.msg)
end
