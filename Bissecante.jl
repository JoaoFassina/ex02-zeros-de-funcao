function Bissecante(enl :: EquacaoNL{T};
                   atol = √eps(T),
                   rtol = √eps(T),
                   max_iter = 10_000,
                   max_time = 30.0,
                   max_eval = 10_000,
                  ) where T <: AbstractFloat


  #aqui vamos usar a bisseccao para diminuir o intervalo, chegando em um ponto inicial mais próximo da raiz
  #a bisseccao tambem exclui as raizes duplas das funcoes, o que "filtra" melhor nosso objetivo
  c=0
  start_time = time()
  f(x) = fun_val(enl, x)
  a = enl.x₀
  fa = f(a)
  ϵ = atol + rtol * abs(fa)
  if abs(fa) ≤ ϵ
    return a, fa, :resolvido, 0, time() - start_time
  end
  b = a + 1
  fb = f(b)
  if abs(fb) ≤ ϵ
    return b, fb, :resolvido, 0, time() - start_time
  end

  while fa * fb ≥ 0
    c += 1
    δ = (b - 0.4*a)
    if δ > 1e5
          return a, fa, :falha, 0, time() - start_time
    end
    if abs(fa) < abs(fb)
          a -= δ
        if c%3 != 0
              fa = f(a)
        end
    else
          b += δ
          if c%3 != 0
              fb = f(b)
        end
    end
  end

  x = (a + b) / 2
  fx = f(x)
  resolvido = abs(fx) ≤ ϵ
  Δt = time() - start_time
  iter = 0
  excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
  proximo = abs(fx) ≤ 0.01 #coloquei esse valor mas não sei se é melhor colocar mais um zero, testa pra ver
  while !(resolvido || excedido || proximo)
    if fa * fx < 0
      b, fb = x, fx
    else
      a, fa = x, fx
    end
    x = (a + b) / 2
    fx = f(x)
    resolvido = abs(fx) ≤ ϵ
    iter += 1
    Δt = time() - start_time
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
    proximo = abs(fx) ≤ 0.01
  end
#analisar se resolveu com a bisseccao
    status = :desconhecido
    if resolvido
      status = :resolvido
    elseif excedido
      if Δt > max_time
      status = :max_time
      else
      status = :max_iter
      end
    fa = f(a)
    fb = f(b)
    elseif proximo
      status = :desconhecido
      while !(resolvido || excedido)
      gx = ((fb - fa) / (b - a))
    	if abs(gx) < atol
      status = :falha
      break
    end
    b, fb = a, fa
    a -= fa / gx
    fa = f(a)
    resolvido = abs(fa) ≤ ϵ
    Δt = time() - start_time
    iter += 1
    excedido = Δt > max_time || iter ≥ max_iter || contador(enl) ≥ max_eval
      end

      if resolvido
        status = :resolvido
      elseif excedido
        if Δt > max_time
          status = :max_time
        else
          status = :max_iter
        end
      end

    end
  return x, fx, status, iter, Δt
end
