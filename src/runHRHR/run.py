from utils.functions import RunModel


model = RunModel()

dose11 = [1]*5
dose12 = [1]*5
dose21 = [1]*5
dose22 = [1]*5

r1 = 10**(-5)
r2 = 10**(-5)

out = model.master_loop_one_tactic(dose11, dose12, dose21, dose22, res_prop_1=r1, res_prop_2=r2)
print(out)