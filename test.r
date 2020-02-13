try(dyn.unload("test"))
compile("test.cpp")
dyn.load("test")

Data = list(i=1)
Parameters = list(p= 0)

# Make object
Obj = MakeADFun(data=Data,
                parameters=Parameters,
                hessian=TRUE,
                DLL="test")

out <- nlminb(Obj$par,Obj$fn,Obj$gr)

rep <- Obj$report()
