library(shiny)
library(visNetwork)


trim_ws <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}

parse_node_set <- function(x) {
  x <- trim_ws(x)
  if (x == "") return(character(0))
  parts <- unlist(strsplit(x, ","))
  parts <- trim_ws(parts)
  parts[parts != ""]
}

parse_edges <- function(text) {
  lines <- unlist(strsplit(text, "\n"))
  lines <- trim_ws(lines)
  lines <- lines[lines != ""]
  
  result <- data.frame(
    edge_id = integer(),
    tail = character(),
    head = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], "->")[[1]]
    
    if (length(parts) != 2) {
      stop(paste("Line", i, "must contain exactly one '->'"))
    }
    
    tail_nodes <- parse_node_set(parts[1])
    head_nodes <- parse_node_set(parts[2])
    
    if (length(tail_nodes) == 0 || length(head_nodes) == 0) {
      stop(paste("Line", i, "must have non-empty tail and head"))
    }
    
    result <- rbind(
      result,
      data.frame(
        edge_id = i,
        tail = paste(tail_nodes, collapse = ","),
        head = paste(head_nodes, collapse = ","),
        stringsAsFactors = FALSE
      )
    )
  }
  
  result
}

get_all_nodes <- function(edge_df) {
  tail_nodes <- unlist(strsplit(paste(edge_df$tail, collapse = ","), ","))
  head_nodes <- unlist(strsplit(paste(edge_df$head, collapse = ","), ","))
  nodes <- unique(c(trim_ws(tail_nodes), trim_ws(head_nodes)))
  nodes[nodes != ""]
}

compute_node_degrees <- function(edge_df) {
  nodes <- get_all_nodes(edge_df)
  
  out_deg <- setNames(rep(0, length(nodes)), nodes)
  in_deg  <- setNames(rep(0, length(nodes)), nodes)
  
  for (i in seq_len(nrow(edge_df))) {
    tail_nodes <- parse_node_set(edge_df$tail[i])
    head_nodes <- parse_node_set(edge_df$head[i])
    
    for (v in tail_nodes) {
      out_deg[v] <- out_deg[v] + 1
    }
    
    for (u in head_nodes) {
      in_deg[u] <- in_deg[u] + 1
    }
  }
  
  data.frame(
    node = nodes,
    out_degree = as.integer(out_deg[nodes]),
    in_degree = as.integer(in_deg[nodes]),
    stringsAsFactors = FALSE
  )
}

compute_node_neighbors <- function(edge_df) {
  nodes <- get_all_nodes(edge_df)
  
  out_neighbors <- setNames(vector("list", length(nodes)), nodes)
  in_neighbors  <- setNames(vector("list", length(nodes)), nodes)
  
  for (v in nodes) {
    out_set <- character(0)
    in_set <- character(0)
    
    for (i in seq_len(nrow(edge_df))) {
      tail_nodes <- parse_node_set(edge_df$tail[i])
      head_nodes <- parse_node_set(edge_df$head[i])
      
      if (v %in% tail_nodes) {
        out_set <- union(out_set, head_nodes)
      }
      
      if (v %in% head_nodes) {
        in_set <- union(in_set, tail_nodes)
      }
    }
    
    out_neighbors[[v]] <- setdiff(out_set, v)
    in_neighbors[[v]]  <- setdiff(in_set, v)
  }
  
  data.frame(
    node = nodes,
    out_neighbors = sapply(nodes, function(v) paste(out_neighbors[[v]], collapse = ",")),
    n_out = sapply(nodes, function(v) length(out_neighbors[[v]])),
    in_neighbors = sapply(nodes, function(v) paste(in_neighbors[[v]], collapse = ",")),
    n_in = sapply(nodes, function(v) length(in_neighbors[[v]])),
    stringsAsFactors = FALSE
  )
}

compute_edge_flow <- function(edge_df) {
  n <- nrow(edge_df)
  
  result <- data.frame(
    edge_id = integer(),
    tail = character(),
    head = character(),
    n_e_in = integer(),
    n_e_out = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n)) {
    Ti <- parse_node_set(edge_df$tail[i])
    Hi <- parse_node_set(edge_df$head[i])
    
    n_in <- 0
    n_out <- 0
    
    for (j in seq_len(n)) {
      if (i == j) next
      
      Tj <- parse_node_set(edge_df$tail[j])
      Hj <- parse_node_set(edge_df$head[j])
      
      if (length(intersect(Hj, Ti)) > 0) {
        n_in <- n_in + 1
      }
      
      if (length(intersect(Tj, Hi)) > 0) {
        n_out <- n_out + 1
      }
    }
    
    result <- rbind(
      result,
      data.frame(
        edge_id = i,
        tail = edge_df$tail[i],
        head = edge_df$head[i],
        n_e_in = n_in,
        n_e_out = n_out,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result
}

compute_dhlrc <- function(edge_df, node_deg_df, edge_flow_df) {
  result <- data.frame(
    edge_id = integer(),
    tail = character(),
    head = character(),
    tail_size = integer(),
    head_size = integer(),
    first_term = numeric(),
    C_e = numeric(),
    M_e = numeric(),
    m_e = numeric(),
    dHLRC = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(edge_df))) {
    T_nodes <- parse_node_set(edge_df$tail[i])
    H_nodes <- parse_node_set(edge_df$head[i])
    
    tail_vals <- node_deg_df$out_degree[match(T_nodes, node_deg_df$node)]
    head_vals <- node_deg_df$in_degree[match(H_nodes, node_deg_df$node)]
    
    tail_safe <- pmax(tail_vals, 1)
    head_safe <- pmax(head_vals, 1)
    
    first_term <- sum(1 / tail_safe) + sum(1 / head_safe)
    
    n_e_in <- edge_flow_df$n_e_in[i]
    n_e_out <- edge_flow_df$n_e_out[i]
    tail_size <- length(T_nodes)
    head_size <- length(H_nodes)
    
    C_e <- n_e_in + n_e_out + (tail_size + head_size) / 2 - 1
    
    all_vals <- c(tail_safe, head_safe)
    M_e <- max(all_vals)
    m_e <- min(all_vals)
    
    raw_score <- first_term + C_e / M_e + C_e / m_e - 1
    
    dhlrc_val <- raw_score
    
    result <- rbind(
      result,
      data.frame(
        edge_id = i,
        tail = edge_df$tail[i],
        head = edge_df$head[i],
        tail_size = tail_size,
        head_size = head_size,
        first_term = first_term,
        C_e = C_e,
        M_e = M_e,
        m_e = m_e,
        dHLRC = dhlrc_val,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result
}

# --- 完美适配无界曲率的可视化逻辑 ---
build_vis_data <- function(edge_df, dhlrc_df) {
  all_nodes <- get_all_nodes(edge_df)
  node_ids <- setNames(seq_along(all_nodes), all_nodes)
  
  nodes_df <- data.frame(
    id = unname(node_ids),
    label = names(node_ids),
    shape = "dot",
    size = 18,
    stringsAsFactors = FALSE
  )
  
  edge_rows <- list()
  k <- 1
  
  for (i in seq_len(nrow(edge_df))) {
    T_nodes <- parse_node_set(edge_df$tail[i])
    H_nodes <- parse_node_set(edge_df$head[i])
    score <- dhlrc_df$dHLRC[dhlrc_df$edge_id == i]
    
    edge_color <- if (score >= 1) {
      "#d73027"  
    } else if (score >= 0) {
      "#fc8d59" 
    } else if (score >= -2) {
      "#91bfdb"  
    } else {
      "#4575b4"  
    }
    
    edge_width <- min(10, 1 + 2 * abs(score))
    
    for (t in T_nodes) {
      for (h in H_nodes) {
        edge_rows[[k]] <- data.frame(
          from = node_ids[[t]],
          to = node_ids[[h]],
          arrows = "to",
          label = paste0("e", i, ": ", round(score, 3)),
          color = edge_color,
          width = edge_width,
          smooth = FALSE,
          stringsAsFactors = FALSE
        )
        k <- k + 1
      }
    }
  }
  
  edges_df <- do.call(rbind, edge_rows)
  list(nodes = nodes_df, edges = edges_df)
}

# -----------------------------
# UI
# -----------------------------

ui <- fluidPage(
  titlePanel("Directed HLRC Explorer & Visualizer (Authentic Unbounded Version)"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "edge_text",
        "Directed hyperedges (Tail -> Head)",
        rows = 10,
        value = paste(
          "A -> X1",
          "A -> X2",
          "A -> X3",
          "Y1 -> B",
          "Y2 -> B",
          "Y3 -> B",
          "A -> B",
          sep = "\n"
        )
      ),
      actionButton("run_calc", "Compute & Visualize"),
      br(), br(),
      tags$p("Edge color legend (Unbounded):"),
      tags$div(style = "color:#d73027; font-weight:bold;", "High Curvature (>= 1) - Strong Community"),
      tags$div(style = "color:#fc8d59; font-weight:bold;", "Moderate Curvature (0 to 1)"),
      tags$div(style = "color:#91bfdb; font-weight:bold;", "Negative Curvature (-2 to 0) - Bridge"),
      tags$div(style = "color:#4575b4; font-weight:bold;", "Deep Negative (< -2) - Severe Bottleneck")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Network Visualization", br(), visNetworkOutput("network_plot", height = "650px")),
        tabPanel("Directed HLRC Results", br(), tableOutput("dhlrc_table")),
        tabPanel("Parsed Hyperedges", br(), tableOutput("edge_table")),
        tabPanel("Node Degrees", br(), tableOutput("degree_table")),
        tabPanel("Edge Flow Neighbors", br(), tableOutput("edge_flow_table"))
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------

server <- function(input, output, session) {
  
  parsed <- eventReactive(input$run_calc, {
    parse_edges(input$edge_text)
  })
  
  node_deg <- eventReactive(input$run_calc, {
    compute_node_degrees(parsed())
  })
  
  node_nb <- eventReactive(input$run_calc, {
    compute_node_neighbors(parsed())
  })
  
  edge_flow <- eventReactive(input$run_calc, {
    compute_edge_flow(parsed())
  })
  
  dhlrc_res <- eventReactive(input$run_calc, {
    # 确保传入的是 node_deg()
    compute_dhlrc(parsed(), node_deg(), edge_flow())
  })
  
  vis_data <- eventReactive(input$run_calc, {
    build_vis_data(parsed(), dhlrc_res())
  })
  
  output$edge_table <- renderTable({
    parsed()
  })
  
  output$degree_table <- renderTable({
    node_deg()
  })
  
  output$neighbor_table <- renderTable({
    node_nb()
  })
  
  output$edge_flow_table <- renderTable({
    edge_flow()
  })
  
  output$dhlrc_table <- renderTable({
    dhlrc_res()
  }, digits = 4)
  
  output$network_plot <- renderVisNetwork({
    vd <- vis_data()
    visNetwork(vd$nodes, vd$edges) %>%
      visEdges(font = list(align = "top")) %>%
      visNodes(font = list(size = 24)) %>%
      visPhysics(stabilization = TRUE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
  })
}

shinyApp(ui = ui, server = server)