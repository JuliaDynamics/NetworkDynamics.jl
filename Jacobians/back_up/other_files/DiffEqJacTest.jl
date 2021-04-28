using DifferentialEquations

function rober(du,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),[0.04,3e7,1e4])
sol = solve(prob, Rodas5())
plot(sol,tspan=(1e-2,1e5),xscale=:log10)


function rober_jac(J,u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  J[1,1] = k₁ * -1
  J[2,1] = k₁
  J[3,1] = 0
  J[1,2] = y₃ * k₃
  J[2,2] = y₂ * k₂ * -2 + y₃ * k₃ * -1
  J[3,2] = y₂ * 2 * k₂
  J[1,3] = k₃ * y₂
  J[2,3] = k₃ * y₂ * -1
  J[3,3] = 0
  nothing
end
f = ODEFunction(rober, jac=rober_jac)
prob_jac = ODEProblem(f,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4))

soljac = solve(prob_jac,Rodas5() )
plot(soljac,tspan=(1e-2,1e5),xscale=:log10)
