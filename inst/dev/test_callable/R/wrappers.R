callable_rpg_hybrid <- function(n, h, z)
    .Call("C_rpg_hybrid", as.double(h), as.double(z), as.integer(n))

callable_rpg_hybrid_fill <- function(h, z)
    .Call("C_rpg_hybrid_fill", as.double(h), as.double(z))

callable_rpg_gamma <- function(n, h, z, trunc=200L)
    replicate(n, .Call("C_rpg_gamma", as.double(h), as.double(z), as.integer(trunc)))

callable_rpg_devroye <- function(n, h, z)
    replicate(n, .Call("C_rpg_devroye", as.integer(h), as.double(z)))

callable_rpg_sp <- function(n, h, z)
    replicate(n, .Call("C_rpg_sp", as.double(h), as.double(z)))
