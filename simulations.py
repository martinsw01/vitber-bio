from montecarlo import monte_carlo
from visualization import show_monomers
from grid_creator import create_grid
import time

#Oppgave 1f)
def simulate_monomers(N, M, N_s, T):
    grid = create_grid(N, M)
    grid, epsilon = monte_carlo(grid, N_s, M, T)
    show_monomers(grid, M)

if __name__ == "__main__":


    t0=time.time()
    simulate_monomers(15, 25, 50000, 500)
    print(time.time()-t0)