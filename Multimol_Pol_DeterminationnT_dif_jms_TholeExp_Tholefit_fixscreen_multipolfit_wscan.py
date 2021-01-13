m=minuit.Minuit(etatot,  screen=2.8990, fix_screen=True)
m.scan()
m.printMode = 1
m.migrad()
m=minuit.Minuit(etatot,  screen=2.8990, fix_screen=True)
m.scan()
m.printMode = 1
m.migrad()
m=minuit.Minuit(etatot,  screen=2.8990, fix_screen=True)
m.scan()
m.printMode = 1
m.migrad()
m=minuit.Minuit(etatot, limit_C=(0.0100,20.0),H=0.001,fix_H=True,limit_S=(0.0100,20.0),limit_single=(0.0100,20.0),limit_pi=(0.0100,20.0), screen=2.8990, fix_screen=True)
m.scan(("C",5,0.0100,20.0),("S",5,0.0100,20.0),("single",5,0.0100,20.0),("pi",5,0.0100,20.0),)
m.printMode = 1
m.migrad()
m=minuit.Minuit(etatot, limit_C=(0.0100,20.0),H=0.001,fix_H=True,limit_S=(0.0100,20.0),limit_single=(0.0100,20.0),limit_pi=(0.0100,20.0), screen=2.8990, fix_screen=True)
m.scan(("C",5,0.0100,20.0),("S",5,0.0100,20.0),("single",5,0.0100,20.0),("pi",5,0.0100,20.0),)
m.printMode = 1
m.migrad()
