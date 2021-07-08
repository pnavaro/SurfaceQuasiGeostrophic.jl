struct TimeSolver

    y0 :: Array{ComplexF64, 2}
    dy :: Array{ComplexF64, 2}

    function TimeSolver( nx, ny)

        y0 = zeros(ComplexF64, nx÷2+1, ny)
        dy = zeros(ComplexF64, nx÷2+1, ny)

        new( y0, dy )

    end

end
