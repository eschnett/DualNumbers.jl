using DualNumbers
using StaticArrays
using Test

@testset "Ring" begin
    T = Int

    randint() = T(rand(-100:100))

    for iter in 1:100
        n = Dual(T(0))
        e = Dual(T(1))
        x = Dual(randint(), randint())
        y = Dual(randint(), randint())
        z = Dual(randint(), randint())
        a = randint()
        b = randint()
        i = rand(1:10)

        @test Dual(a).primal === a
        @test Dual(a).dual === 0
        @test Dual(a, b).primal === a
        @test Dual(a, b).dual === b

        @test +x === x
        @test -(-x) === x

        @test n + x === x
        @test x + n === x

        @test -x === n - x
        @test x - x === n

        @test x + y === y + x

        @test (x + y) + z === x + (y + z)

        @test n * x === n
        @test x * n === n

        @test e * x === x
        @test x * e === x

        @test x * y === y * x
        @test (x * y) * z === x * (y * z)

        @test x^0 == e
        @test x^1 == x
        @test x^i == x * x^(i - 1)

        @test x + a === a + x

        @test Dual(a) == a
        @test (Dual(a) == b) == (a == b)
        @test (Dual(a) < b) == (a < b)
        @test (Dual(a) > b) == (a > b)
    end
end

@testset "Field" begin
    T = Rational{Int128}

    randint() = T(rand(-100:100))
    randintnz() =
        while true
            i = randint()
            i != 0 && return i
        end
    randrat() = randint() // randintnz()

    for iter in 1:100
        n = Dual(T(0))
        e = Dual(T(1))
        x = Dual(randrat(), randrat())
        y = Dual(randrat(), randrat())
        z = Dual(randrat(), randrat())
        a = randrat()
        b = randrat()

        if x != 0
            @test inv(x) === 1 / x
            @test inv(inv(x)) === x

            @test e / x === inv(x)
            @test x \ y === y / x

            if y != 0
                @test x / y === inv(y / x)
                @test x \ y === inv(y \ x)
            end
        end
    end
end

@testset "Real" begin
    T = Float64

    randfloat() = rand(T) * 200 - 100

    for iter in 1:100
        n = Dual(T(0))
        e = Dual(T(1))
        x = Dual(randfloat(), randfloat())
        y = Dual(randfloat(), randfloat())
        z = Dual(randfloat(), randfloat())
        a = randfloat()
        b = randfloat()

        @test (abs(x)^a).primal ≈ abs(x).primal^a
        @test cos(x).primal ≈ cos(x.primal)
        @test sin(x).primal ≈ sin(x.primal)
        @test sqrt(abs(x)).primal ≈ sqrt(abs(x).primal)
    end
end

@testset "Basic Derivatives" begin
    T = Float64

    randfloat() = rand(T)

    for iter in 1:100
        x = randfloat()
        a = randfloat()

        @test derivative(x -> x + a, x) ≈ 1
        @test derivative(x -> a * x, x) ≈ a
        @test derivative(x -> x^2, x) ≈ 2x
        @test derivative(cos, x) ≈ -sin(x)
        @test derivative(sin, x) ≈ cos(x)
        @test derivative(sqrt, x) ≈ 1 / 2sqrt(x)
    end
end

@testset "Derivative laws" begin
    T = Float64

    randfloat() = rand(T)

    for iter in 1:100
        x = randfloat()

        c0 = randfloat()
        c1 = randfloat()
        fs = [x -> x + c0, x -> c1 * x, x -> x^2, cos, sin, x -> sqrt(x + 1)]
        i = identity
        f = rand(fs)
        g = rand(fs)
        h = rand(fs)

        a = randfloat()

        @test (g ∘ f)(x) ≈ g(f(x))

        @test derivative(x -> f(x) + g(x), x) ≈ derivative(f, x) + derivative(g, x)
        @test derivative(x -> a * f(x), x) ≈ a * derivative(f, x)

        @test derivative(identity ∘ f, x) ≈ derivative(f, x)
        @test derivative(f ∘ identity, x) ≈ derivative(f, x)

        @test derivative(identity, x) ≈ 1
        f′(x) = derivative(f, x)
        g′(x) = derivative(g, x)
        @test derivative(g ∘ f, x) ≈ g′(f(x)) * f′(x)
    end
end

@testset "Derivatives of multi-valued functions" begin
    T = Float64

    randfloat() = rand(T)

    for iter in 1:100
        xs = SVector(randfloat(), randfloat(), randfloat())

        c0 = randfloat()
        c1 = randfloat()
        fs = [x -> x + c0, x -> c1 * x, x -> x^2, cos, sin, x -> sqrt(x + 1)]
        f1 = rand(fs)
        f2 = rand(fs)
        f3 = rand(fs)
        f(xs) = f1(xs[1]) + f2(xs[2]) * f3(xs[3])
        g(xs) = f1(xs[1]) * f2(xs[2]) * f3(xs[3])

        a = randfloat()

        @test derivative(f, xs, 1) ≈ derivative(f1, xs[1])
        @test derivative(f, xs, 2) ≈ derivative(f2, xs[2]) * f3(xs[3])
        @test derivative(f, xs, 3) ≈ f2(xs[2]) * derivative(f3, xs[3])

        for i in 1:3
            @test derivative(xs -> f(xs) + g(xs), xs, i) ≈ derivative(f, xs, i) + derivative(g, xs, i)
            @test derivative(xs -> a * f(xs), xs, i) ≈ a * derivative(f, xs, i)
        end
    end
end
