plot_ini_ens <- function(ini_ens, reference, observations, 
                         print_thickness_plot = TRUE, print_bed_plot = TRUE, 
                         print_alpha_plot = FALSE, print_friction_plot = TRUE) {
  
  domain <- reference$domain
  
  ############ ggplot for initial ensemble ###########
  ## ggplot for ice thickness ##
  
  par(mfrow = c(2,2))
  
  if (print_thickness_plot) {
  ini_thickness <- data.frame(ini_ens[1:J, ], mean = rowMeans(ini_ens[1:J, ]),
                              domain = domain/1000)
  reference_thickness <- data.frame(reference$all_thicknesses, domain = domain/1000)
  
  # Grounding line
  ini_GL_ref <- reference$grounding_line[1]
  ini_GL <- gl_migrate(H = rowMeans(ini_ens[1:J, ]), b = rowMeans(ini_ens[(J+1):(2*J), ]))
  ini_GL <- ini_GL / J * domain[length(domain)] / 1000
  
  ini_thickness_plot <- ggplot() + 
    # geom_ribbon(data = filtered_thickness, 
    #                    aes(domain, filtered_thickness[, t], ymin = thickness_low, ymax = thickness_up),
    #                    fill = "red", alpha = 0.2) 
    theme_bw() +
    geom_line(data = ini_thickness, 
              aes(domain, ini_thickness[, 1]), colour = "lightblue") +
    geom_line(data = ini_thickness, 
              aes(domain, ini_thickness[, 2]), colour = "pink") +
    geom_line(data = ini_thickness, 
              aes(domain, ini_thickness[, 3]), colour = "plum") +
    geom_line(data = reference_thickness, 
              aes(domain, reference_thickness[, 1]), colour = "black") +
    geom_line(data = ini_thickness,
              aes(domain, mean), colour = "dodgerblue") +
    geom_vline(xintercept = ini_GL_ref, linetype="dotted", 
               colour = "black", size = 0.5) +
    geom_vline(xintercept = ini_GL, linetype="dashed",
               colour = "dodgerblue", size = 0.5) +
    coord_cartesian(xlim = c(0, 600)) +
    xlab("x (km)") + ylab("Ice thickness (m)") 
  
  
    print(ini_thickness_plot)
  }
  
  
  ## Initial bedrock
  if (print_bed_plot) {
    dx <- domain[2] - domain[1] ## Fix this later so that dx is an output from model run
    bed_obs <- data.frame(bed = observations$bed_obs, 
                          locations = observations$bed_obs_locations * dx / 1000)
    
    ini_bed <- data.frame(ini_ens[(J+1):(2*J), ], mean = rowMeans(ini_ens[(J+1):(2*J), ]),
                          domain = domain/1000)
    reference_bed <- data.frame(bed = reference$bedrock, domain = domain/1000)
    
    ini_bed_plot <- ggplot() + 
      theme_bw() +
      geom_line(data = ini_bed, 
                aes(domain, ini_bed[, 1]), colour = "lightblue") +
      geom_line(data = ini_bed, 
                aes(domain, ini_bed[, 2]), colour = "pink") +
      geom_line(data = ini_bed, 
                aes(domain, ini_bed[, 3]), colour = "plum") +
      geom_line(data = reference_bed, 
                aes(domain, bed), colour = "black") +
      geom_line(data = ini_bed,
                aes(domain, mean), colour = "dodgerblue") +
      geom_point(data = bed_obs, 
                 aes(locations, bed), colour = "cyan") +
      geom_vline(xintercept = ini_GL_ref, linetype="dotted", 
                 colour = "black", size = 0.5) +
      geom_vline(xintercept = ini_GL, linetype="dashed",
                 colour = "dodgerblue", size = 0.5) +
      coord_cartesian(xlim = c(0, 600)) +
      xlab("x (km)") + ylab("Bed elevation (m)") 
  
    print(ini_bed_plot)
  }
  
  ## Initial Alpha ##
  if (print_alpha_plot) {
    secpera <- 31556926 #seconds per annum
    fric_scale <- secpera^(-1/3) / 1e6
    ini_alpha <- data.frame(ini_ens[(2*J+1):(3*J), ] * fric_scale, 
                            mean = rowMeans(ini_ens[(2*J+1):(3*J), ]) * fric_scale,
                            domain = domain/1000)
    reference_alpha <- data.frame(alpha = log10(reference$friction_coef) * fric_scale, 
                                  domain = domain/1000)
    
    ini_alpha_plot <- ggplot() + 
      theme_bw() +
      geom_line(data = ini_alpha, 
                aes(domain, ini_alpha[, 1]), colour = "lightblue") +
      geom_line(data = ini_alpha, 
                aes(domain, ini_alpha[, 2]), colour = "pink") +
      geom_line(data = ini_alpha, 
                aes(domain, ini_alpha[, 3]), colour = "plum") +
      geom_line(data = reference_alpha, 
                aes(domain, alpha), colour = "black") +
      geom_line(data = ini_alpha,
                aes(domain, mean), colour = "dodgerblue", linewidth = 0.75) +
      geom_vline(xintercept = ini_GL_ref, linetype="dotted", 
                 colour = "black", linewidth = 0.75) +
      geom_vline(xintercept = ini_GL, linetype="dashed",
                 colour = "dodgerblue", linewidth = 0.75) +
      coord_cartesian(xlim = c(230, 370), ylim = c(1.75e-08, 2.25e-08)) +
      xlab("x (km)") + ylab(expression(alpha~(M~Pa~m^-{1/3}~a^{1/3}))) 
    
    print(ini_alpha_plot)
  }
  
  
  ## Initial basal friction coef ##
  if (print_friction_plot) {
    secpera <- 31556926 #seconds per annum
    fric_scale <- secpera^(-1/3) / 1e6
    
    ini_friction <- data.frame(10^ini_ens[(2*J+1):(3*J), ] * fric_scale, 
                               mean = rowMeans(10^ini_ens[(2*J+1):(3*J), ]) * fric_scale,
                               domain = domain/1000)
    reference_friction <- data.frame(friction = reference$friction_coef * fric_scale, 
                                     domain = domain/1000)
    
    ini_friction_plot <- ggplot() + 
      theme_bw() +
      geom_line(data = ini_friction, 
                aes(domain, ini_friction[, 1]), colour = "lightblue") +
      geom_line(data = ini_friction, 
                aes(domain, ini_friction[, 2]), colour = "pink") +
      geom_line(data = ini_friction, 
                aes(domain, ini_friction[, 3]), colour = "plum") +
      geom_line(data = reference_friction, 
                aes(domain, friction), colour = "black") +
      geom_line(data = ini_friction,
                aes(domain, mean), colour = "dodgerblue", linewidth = 0.75) +
      geom_vline(xintercept = ini_GL_ref, linetype="dotted", 
                 colour = "black", linewidth = 0.75) +
      geom_vline(xintercept = ini_GL, linetype="dashed",
                 colour = "dodgerblue", linewidth = 0.75) +
      coord_cartesian(xlim = c(230, 370), ylim = c(0, 0.05)) +
      xlab("x (km)") + ylab(expression(C~(M~Pa~m^-{1/3}~a^{1/3}))) 
  
    print(ini_friction_plot)
  }
  
}