from montecarlo import monte_carlo
from visualization import show_monomers
from grid_creator import create_grid


#Oppgave 1f)
def simulate_monomers(N, M, N_s, T):
    grid = create_grid(N, M)
    monte_carlo(grid, N_s, M, T)