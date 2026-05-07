import tkinter as tk
from file_selection import root



class App(tk.Tk) : 
    def __init__(self, *args, **kwargs):

        tk.TK.__init__(self, *args, **kwargs)

        container = tk.Frame(self)
        container.pack(site="top", fill = "both", expand = True)

        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (plink_select, sim_covariates, sim_phenotypes) :
            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row = 0, column = 0, sticky = "nsew")


# ask user for number of individuals
num_individ_lab = tk.Label(root, text="Enter number of individuals:")
num_individ = tk.Entry(root)
# ask user for number of SNPs
num_snps_lab = tk.Label(root, text= "Enter number of SNPs:")
num_snps = tk.Entry(root)

# ask user for number of populations
num_pop_lab = tk.Label(root, text="Enter number of populations:")
num_pops = tk.Entry(root)

# ask user for number of sites
num_sites_lab = tk.Label(root, text="Enter number of sites:")
num_sites = tk.Entry(root)

display = tk.Label(root, text="make a selection and click update") 

# a functino that updates the display string with the user selected site_effects and pop_dist
def update_string():
    display.config(text="You selected " + pop_dist.get() + " and " + site_effects.get())

# ask user select distribution of populations from dropdown menu
pop_dist_lab = tk.Label(root, text="Select distribution of populations:")
pop_dist = tk.StringVar(root)
pop_dist.set("Select distribution")
pop_dist_menu = tk.OptionMenu(root, pop_dist, "Equal", "IID", "Het")



site_effects_lab = tk.Label(root, text = "Choose site effects:")
# ask user to select type of site effects from dropdown menu
site_effects = tk.StringVar(root)
site_effects.set("Select site effects")
# add a lambda function to get site_effects
site_effects_menu = tk.OptionMenu(root, site_effects, "None", "Random", "Fixed")


update_button = tk.Button(root, text="Update", command=update_string)




# pack elements in the order that they were defined
plink_lab.pack()
plink_filepath.pack()

num_individ_lab.pack()
num_individ.pack()

num_snps_lab.pack()
num_snps.pack()

num_pop_lab.pack()
num_pops.pack()

num_sites_lab.pack()
num_sites.pack()

pop_dist_lab.pack()
pop_dist_menu.pack()

site_effects_lab.pack()
site_effects_menu.pack()

display.pack()
update_button.pack()

root.mainloop()
