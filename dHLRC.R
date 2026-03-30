library(shiny)
library(visNetwork)
library(DT)

# ==========================================
# 1. Core Functions
# ==========================================

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
    parts <- strsplit(lines[i], "->", fixed = TRUE)[[1]]
    
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

# ==========================================
# 2. Degree-based quantities
# ==========================================

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

# ==========================================
# 3. Neighborhood-based quantities
# ==========================================

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

# ==========================================
# 4. Edge flow counts
# Exclude self-edge to avoid self double counting
# ==========================================

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

# ==========================================
# 5. Main dHLRC calculator
# Simultaneously compute both original and normalized scores
# ==========================================

compute_dhlrc <- function(edge_df, node_deg_df, node_nb_df, edge_flow_df,
                          metric_type = "neighbor") {
  
  result <- data.frame(
    edge_id = integer(),
    tail = character(),
    head = character(),
    dHLRC_original = numeric(),
    dHLRC_normalized = numeric(),
    gap = numeric(),
    tail_size = integer(),
    head_size = integer(),
    metric_type = character(),
    first_term = numeric(),
    C_original = numeric(),
    C_normalized = numeric(),
    M_e = numeric(),
    m_e = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(edge_df))) {
    T_nodes <- parse_node_set(edge_df$tail[i])
    H_nodes <- parse_node_set(edge_df$head[i])
    
    if (metric_type == "degree") {
      tail_vals <- node_deg_df$out_degree[match(T_nodes, node_deg_df$node)]
      head_vals <- node_deg_df$in_degree[match(H_nodes, node_deg_df$node)]
    } else {
      tail_vals <- node_nb_df$n_out[match(T_nodes, node_nb_df$node)]
      head_vals <- node_nb_df$n_in[match(H_nodes, node_nb_df$node)]
    }
    
    tail_safe <- pmax(tail_vals, 1)
    head_safe <- pmax(head_vals, 1)
    
    first_term <- 0.5 * (sum(1 / tail_safe) + sum(1 / head_safe))
    
    n_e_in <- edge_flow_df$n_e_in[i]
    n_e_out <- edge_flow_df$n_e_out[i]
    tail_size <- length(T_nodes)
    head_size <- length(H_nodes)
    
    C_original <- n_e_in + n_e_out + (tail_size + head_size) / 2 - 1
    C_normalized <- (n_e_in + n_e_out) / 2 + (tail_size + head_size) / 4 - 1
    
    all_vals <- c(tail_safe, head_safe)
    M_e <- max(all_vals)
    m_e <- min(all_vals)
    
    dhlrc_original <- first_term + C_original / M_e + C_original / m_e - 1
    dhlrc_normalized <- first_term + C_normalized / M_e + C_normalized / m_e - 1
    
    result <- rbind(
      result,
      data.frame(
        edge_id = i,
        tail = edge_df$tail[i],
        head = edge_df$head[i],
        dHLRC_original = dhlrc_original,
        dHLRC_normalized = dhlrc_normalized,
        gap = dhlrc_original - dhlrc_normalized,
        tail_size = tail_size,
        head_size = head_size,
        metric_type = metric_type,
        first_term = first_term,
        C_original = C_original,
        C_normalized = C_normalized,
        M_e = M_e,
        m_e = m_e,
        stringsAsFactors = FALSE
      )
    )
  }
  
  result
}

# ==========================================
# 6. Visualization
# score_type: "original" or "normalized"
# ==========================================

build_vis_data <- function(edge_df, dhlrc_df, pos_thresh, neg_thresh,
                           score_type = "original") {
  all_nodes <- get_all_nodes(edge_df)
  node_ids <- setNames(seq_along(all_nodes), all_nodes)
  
  nodes_df <- data.frame(
    id = unname(node_ids),
    label = names(node_ids),
    shape = "dot",
    size = 18,
    stringsAsFactors = FALSE
  )
  
  score_col <- if (score_type == "normalized") {
    "dHLRC_normalized"
  } else {
    "dHLRC_original"
  }
  
  edge_rows <- list()
  k <- 1
  
  for (i in seq_len(nrow(edge_df))) {
    T_nodes <- parse_node_set(edge_df$tail[i])
    H_nodes <- parse_node_set(edge_df$head[i])
    score <- dhlrc_df[[score_col]][dhlrc_df$edge_id == i][1]
    
    edge_color <- if (score >= pos_thresh) {
      "#d73027"
    } else if (score >= 0) {
      "#fc8d59"
    } else if (score >= neg_thresh) {
      "#91bfdb"
    } else {
      "#4575b4"
    }
    
    edge_width <- if (score < 0) {
      min(10, 1 + 3 * abs(score))
    } else {
      min(8, 1 + score)
    }
    
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

# ==========================================
# 7. User Interface
# ==========================================

ui <- fluidPage(
  titlePanel("Directed HLRC Explorer & Visualizer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      textAreaInput(
        "edge_text",
        "Directed hyperedges (Tail -> Head)",
        rows = 10,
        value = paste(
          "A -> X1", "A -> X2", "A -> X3",
          "Y1 -> B", "Y2 -> B", "Y3 -> B",
          "A -> B",
          sep = "\n"
        )
      ),
      
      actionButton("run_calc", "Compute & Visualize", class = "btn-primary"),
      
      hr(),
      tags$h4("Display Settings"),
      
      radioButtons(
        "metric_type",
        "Node-level metric:",
        choices = c(
          "Neighborhood-based (HLRC-style)" = "neighbor",
          "Degree-based" = "degree"
        ),
        selected = "neighbor"
      ),
      
      radioButtons(
        "display_score",
        "Network plot uses:",
        choices = c(
          "Original directed score" = "original",
          "Normalized score" = "normalized"
        ),
        selected = "original"
      ),
      
      hr(),
      tags$h4("Plot Thresholds"),
      
      sliderInput(
        "pos_thresh",
        "Positive threshold (red):",
        min = 0.1, max = 10.0, value = 2.0, step = 0.5
      ),
      sliderInput(
        "neg_thresh",
        "Negative threshold (dark blue):",
        min = -2.0, max = 0.0, value = -0.5, step = 0.1
      ),
      
      br(),
      uiOutput("compact_note"),
      br(),
      uiOutput("plot_legend")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Network Visualization", br(), visNetworkOutput("network_plot", height = "650px")),
        tabPanel("dHLRC Results", br(), DTOutput("dhlrc_table")),
        tabPanel("Parsed Hyperedges", br(), DTOutput("edge_table")),
        tabPanel("Node Degrees", br(), DTOutput("degree_table")),
        tabPanel("Node Neighborhoods", br(), DTOutput("neighbor_table")),
        tabPanel("Edge Flow Neighbors", br(), DTOutput("edge_flow_table"))
      )
    )
  )
)

# ==========================================
# 8. Server Logic
# ==========================================

server <- function(input, output, session) {
  
  parsed <- eventReactive(input$run_calc, {
    parse_edges(input$edge_text)
  }, ignoreNULL = FALSE)
  
  node_deg <- eventReactive(input$run_calc, {
    compute_node_degrees(parsed())
  }, ignoreNULL = FALSE)
  
  node_nb <- eventReactive(input$run_calc, {
    compute_node_neighbors(parsed())
  }, ignoreNULL = FALSE)
  
  edge_flow <- eventReactive(input$run_calc, {
    compute_edge_flow(parsed())
  }, ignoreNULL = FALSE)
  
  dhlrc_res <- eventReactive(input$run_calc, {
    compute_dhlrc(
      edge_df = parsed(),
      node_deg_df = node_deg(),
      node_nb_df = node_nb(),
      edge_flow_df = edge_flow(),
      metric_type = input$metric_type
    )
  }, ignoreNULL = FALSE)
  
  vis_data <- reactive({
    req(dhlrc_res())
    build_vis_data(
      edge_df = parsed(),
      dhlrc_df = dhlrc_res(),
      pos_thresh = input$pos_thresh,
      neg_thresh = input$neg_thresh,
      score_type = input$display_score
    )
  })
  
  output$compact_note <- renderUI({
    metric_txt <- if (input$metric_type == "neighbor") {
      "Node term uses directional neighborhood sizes."
    } else {
      "Node term uses directional degrees."
    }
    
    plot_txt <- if (input$display_score == "normalized") {
      "The network plot is currently colored by the normalized score."
    } else {
      "The network plot is currently colored by the original directed score."
    }
    
    tags$div(
      style = "padding:8px; background:#f7f7f7; border:1px solid #ddd; font-size:13px;",
      tags$p(tags$b("Current setting"), style = "margin-bottom:6px;"),
      tags$p(metric_txt, style = "margin-bottom:4px;"),
      tags$p("The results table always shows both original and normalized scores.", style = "margin-bottom:4px;"),
      tags$p(plot_txt, style = "margin-bottom:0px;")
    )
  })
  
  output$plot_legend <- renderUI({
    tags$div(
      style = "padding:8px; background:#fafafa; border:1px solid #ddd; font-size:13px;",
      tags$p(tags$b("Legend for network plot"), style = "margin-bottom:6px;"),
      tags$div(style = "color:#d73027; font-weight:bold;",
               sprintf(">= %s : strong positive", input$pos_thresh)),
      tags$div(style = "color:#fc8d59; font-weight:bold;",
               sprintf("0 to %s : moderate positive", input$pos_thresh)),
      tags$div(style = "color:#91bfdb; font-weight:bold;",
               sprintf("%s to 0 : mild negative", input$neg_thresh)),
      tags$div(style = "color:#4575b4; font-weight:bold;",
               sprintf("< %s : strong negative", input$neg_thresh))
    )
  })
  
  output$edge_table <- renderDT({
    datatable(
      parsed(),
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })
  
  output$degree_table <- renderDT({
    datatable(
      node_deg(),
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })
  
  output$neighbor_table <- renderDT({
    datatable(
      node_nb(),
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })
  
  output$edge_flow_table <- renderDT({
    datatable(
      edge_flow(),
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE
      )
    )
  })
  
  output$dhlrc_table <- renderDT({
    df <- dhlrc_res()
    
    num_cols <- sapply(df, is.numeric)
    df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
    
    datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(
        scrollX = TRUE,
        pageLength = 10,
        autoWidth = TRUE,
        dom = "tip"
      )
    )
  })
  
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
